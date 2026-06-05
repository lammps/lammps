/*  Style registration generator for LAMMPS
 * (c) 2026 LAMMPS development team <developers@lammps.org>
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * * Redistributions of source code must retain the above copyright
 *   notice, this list of conditions and the following disclaimer.
 * * Redistributions in binary form must reproduce the above copyright
 *   notice, this list of conditions and the following disclaimer in the
 *   documentation and/or other materials provided with the distribution.
 * * Neither the mentioning of the LAMMPS developers nor the names of its
 *   members may be used to endorse or promote products derived from this
 *   software without specific prior written permission.
 *
 *   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 *   AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 *   IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 *   ARE DISCLAIMED. IN NO EVENT SHALL THE LAMMPS DEVELOPERS BE LIABLE FOR ANY
 *   DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 *   (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 *   LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 *   ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 *   (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
 *   THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

/* Alternative to using awk in src/Make.sh.  This keeps the make-based portable
 * (same as fastdep.c, which is also compiled with "cc" during the build).
 *
 * It parses the "XxxStyle(key,Class)" marker lines inside "#ifdef
 * XXX_CLASS" blocks of the LAMMPS style headers and emits the global
 * style-registry registration translation units.
 *
 * Two modes:
 *   stylegen source <macro> <base> <accessor> <regfunc> <onearg> "<includes>" \
 *                   -- file.h ...
 *      Emit a complete style_<x>.cpp to stdout: include list, creator<>
 *      template (one-arg or three-arg signature), and a register function with
 *      one add_builtin() call per parsed marker.  Preprocessor directives
 *      nested inside the marker block (e.g. "#ifdef LMP_KOKKOS_GPU") are copied
 *      through verbatim so build-config-dependent markers stay guarded.
 *
 *   stylegen package <macro> <member> -- file.h ...
 *      Emit "package_styles().<member>[\"key\"] = \"<pkg>\";" lines, where the
 *      package name is the parent directory of each header file.
 */

#include <stddef.h>
#include <stdio.h>
#include <string.h>

#define MAXLINE 4096
#define MAXTOK 512

/* copy the range [b,e) into dst, trimming leading/trailing spaces and tabs */

static void copy_trim(char *dst, const char *b, const char *e)
{
    size_t n;
    while ((b < e) && ((*b == ' ') || (*b == '\t'))) ++b;
    while ((e > b) && ((e[-1] == ' ') || (e[-1] == '\t'))) --e;
    n = (size_t) (e - b);
    if (n >= MAXTOK) n = MAXTOK - 1;
    memcpy(dst, b, n);
    dst[n] = '\0';
}

/* parse "(key , Class)" starting at the '(' character.  Splits the parenthesized
 * text on its first comma; Class may itself contain commas (template arguments).
 * Returns 1 on success and fills key and cls, 0 if malformed. */

static int parse_marker(const char *paren, char *key, char *cls)
{
    const char *close, *comma;
    if (*paren != '(') return 0;
    close = strchr(paren, ')');
    if (!close) return 0;
    comma = memchr(paren + 1, ',', (size_t) (close - (paren + 1)));
    if (!comma) return 0;
    copy_trim(key, paren + 1, comma);
    copy_trim(cls, comma + 1, close);
    return 1;
}

/* return 1 if the line is "^#ifdef[ \t]+<NAME>_CLASS" (NAME uppercase/underscore) */

static int is_class_ifdef(const char *line)
{
    const char *p, *e;
    if (strncmp(line, "#ifdef", 6) != 0) return 0;
    p = line + 6;
    if ((*p != ' ') && (*p != '\t')) return 0;
    while ((*p == ' ') || (*p == '\t')) ++p;
    e = p;
    while (*e && (*e != ' ') && (*e != '\t') && (*e != '\n') && (*e != '\r')) ++e;
    if ((e - p) >= 6 && (strncmp(e - 6, "_CLASS", 6) == 0)) return 1;
    return 0;
}

/* if the line is a preprocessor directive, return a pointer to the directive
 * keyword (after '#' and any whitespace); otherwise return NULL */

static const char *directive_kw(const char *line)
{
    const char *p = line;
    while ((*p == ' ') || (*p == '\t')) ++p;
    if (*p != '#') return NULL;
    ++p;
    while ((*p == ' ') || (*p == '\t')) ++p;
    return p;
}

static int starts_with(const char *s, const char *prefix)
{
    return strncmp(s, prefix, strlen(prefix)) == 0;
}

/* print a directive line with leading whitespace removed and a single newline */

static void emit_directive(const char *line)
{
    const char *p = line;
    size_t n;
    while ((*p == ' ') || (*p == '\t')) ++p;
    n = strlen(p);
    while ((n > 0) && ((p[n - 1] == '\n') || (p[n - 1] == '\r'))) --n;
    printf("%.*s\n", (int) n, p);
}

/* parse one style header and emit add_builtin() calls (source mode).  Runs a
 * small state machine over the "#ifdef XXX_CLASS" block: 0 = before block,
 * 1 = inside block, 2 = after block (the #else class-definition branch). */

static void process_source_file(const char *path, const char *macro, const char *accessor)
{
    FILE *fp;
    char line[MAXLINE];
    int state = 0, depth = 0;
    size_t maclen = strlen(macro);

    fp = fopen(path, "r");
    if (!fp) {
        fprintf(stderr, "stylegen: cannot open %s\n", path);
        return;
    }

    while (fgets(line, sizeof(line), fp)) {
        const char *kw;
        const char *q;

        if (state == 0) {
            if (is_class_ifdef(line)) {
                state = 1;
                depth = 1;
            }
            continue;
        }
        if (state == 2) continue;

        /* state == 1: inside the marker block */
        kw = directive_kw(line);
        if (kw) {
            if (starts_with(kw, "if")) {
                ++depth;
                emit_directive(line);
            } else if (starts_with(kw, "elif")) {
                if (depth > 1) emit_directive(line);
            } else if (starts_with(kw, "else")) {
                if (depth == 1) state = 2;
                else emit_directive(line);
            } else if (starts_with(kw, "endif")) {
                --depth;
                if (depth == 0) state = 2;
                else emit_directive(line);
            }
            /* any other directive (#include, #define, ...) is ignored */
            continue;
        }

        /* not a directive: match "^[ \t]*<macro>(" */
        q = line;
        while ((*q == ' ') || (*q == '\t')) ++q;
        if ((strncmp(q, macro, maclen) == 0) && (q[maclen] == '(')) {
            char key[MAXTOK], cls[MAXTOK];
            if (parse_marker(q + maclen, key, cls))
                printf("  %s().add_builtin(\"%s\", &creator<%s>);\n", accessor, key, cls);
        }
    }
    fclose(fp);
}

/* parse one package style header and emit package_styles() map entries.  Unlike
 * source mode this scans every line for "<macro>(" (the marker only ever appears
 * inside the block, so no state machine is needed) and records key -> package. */

static void process_package_file(const char *path, const char *macro, const char *member,
                                 const char *pkg)
{
    FILE *fp;
    char line[MAXLINE];
    char needle[MAXTOK];
    size_t maclen = strlen(macro);

    if (maclen + 2 > sizeof(needle)) return;
    memcpy(needle, macro, maclen);
    needle[maclen] = '(';
    needle[maclen + 1] = '\0';

    fp = fopen(path, "r");
    if (!fp) {
        fprintf(stderr, "stylegen: cannot open %s\n", path);
        return;
    }

    while (fgets(line, sizeof(line), fp)) {
        char *m = strstr(line, needle);
        if (m) {
            char key[MAXTOK], cls[MAXTOK];
            if (parse_marker(m + maclen, key, cls))
                printf("  package_styles().%s[\"%s\"] = \"%s\";\n", member, key, pkg);
        }
    }
    fclose(fp);
}

/* write the parent directory of path into buf (or "." if there is no slash) */

static const char *dir_name(const char *path, char *buf, size_t bufsz)
{
    const char *slash = strrchr(path, '/');
    size_t n;
    if (!slash) {
        buf[0] = '.';
        buf[1] = '\0';
        return buf;
    }
    n = (size_t) (slash - path);
    if (n >= bufsz) n = bufsz - 1;
    memcpy(buf, path, n);
    buf[n] = '\0';
    return buf;
}

/* print "#include" lines for each whitespace-separated token in s */

static void print_includes(const char *s)
{
    const char *p = s;
    while (*p) {
        const char *start;
        while ((*p == ' ') || (*p == '\t')) ++p;
        if (!*p) break;
        start = p;
        while (*p && (*p != ' ') && (*p != '\t')) ++p;
        printf("#include \"%.*s\"\n", (int) (p - start), start);
    }
}

int main(int argc, char **argv)
{
    int i, sep = -1;
    const char *mode;

    if (argc < 2) {
        fprintf(stderr, "usage: %s source|package ... -- files...\n", argv[0]);
        return 1;
    }
    mode = argv[1];
    for (i = 2; i < argc; ++i) {
        if (strcmp(argv[i], "--") == 0) {
            sep = i;
            break;
        }
    }
    if (sep < 0) {
        fprintf(stderr, "stylegen: missing '--' separator before file list\n");
        return 1;
    }

    if (strcmp(mode, "source") == 0) {
        const char *macro, *base, *accessor, *regfunc, *onearg, *includes;
        if (sep < 8) {
            fprintf(stderr, "stylegen source: not enough arguments\n");
            return 1;
        }
        macro = argv[2];
        base = argv[3];
        accessor = argv[4];
        regfunc = argv[5];
        onearg = argv[6];
        includes = argv[7];

        printf("// auto-generated by the LAMMPS build system from style markers. DO NOT EDIT.\n\n");
        printf("#include \"lammps.h\"\n");
        printf("#include \"creator_registry.h\"\n");
        print_includes(includes);
        printf("\n");
        for (i = sep + 1; i < argc; ++i) printf("#include \"%s\"\n", argv[i]);
        printf("\n");
        printf("namespace LAMMPS_NS {\n\n");
        if (strcmp(onearg, "1") == 0) {
            printf("template <typename T> static %s *creator(LAMMPS *lmp)\n", base);
            printf("{\n  return new T(lmp);\n}\n");
        } else {
            printf("template <typename T> static %s *creator(LAMMPS *lmp, int narg, char **arg)\n", base);
            printf("{\n  return new T(lmp, narg, arg);\n}\n");
        }
        printf("\n");
        printf("void %s()\n{\n", regfunc);
        for (i = sep + 1; i < argc; ++i) process_source_file(argv[i], macro, accessor);
        printf("}\n\n");
        printf("}    // namespace LAMMPS_NS\n");
        return 0;
    }

    if (strcmp(mode, "package") == 0) {
        const char *macro, *member;
        if (sep < 4) {
            fprintf(stderr, "stylegen package: not enough arguments\n");
            return 1;
        }
        macro = argv[2];
        member = argv[3];
        for (i = sep + 1; i < argc; ++i) {
            char pkg[MAXTOK];
            dir_name(argv[i], pkg, sizeof(pkg));
            process_package_file(argv[i], macro, member, pkg);
        }
        return 0;
    }

    fprintf(stderr, "stylegen: unknown mode '%s'\n", mode);
    return 1;
}
