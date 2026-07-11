/*  Fast dependency generator for LAMMPS
 * (c) 2016 Axel Kohlmeyer <akohlmey@gmail.com>
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
 * * Neither the name of the <organization> nor the
 *   names of its contributors may be used to endorse or promote products
 *   derived from this software without specific prior written permission.
 *
 *   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 *   AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 *   IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 *   ARE DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
 *   DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 *   (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 *   LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 *   ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 *   (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
 *   THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>

const char version[] = "2.1";

/* store list of accepted extensions for object targets */
static const char *extensions[] = {".cpp", ".c", ".cu"};
static const int numextensions = sizeof(extensions) / sizeof(const char *);

/* strdup() is not part of ANSI C. provide a replacement for portability */
static char *my_strdup(const char *src)
{
    int len;
    char *ptr;

    if (src == NULL) return NULL;
    len = strlen(src);
    ptr = (char *)malloc(len + 1);
    if (ptr) memcpy(ptr, src, len + 1);
    return ptr;
}

/************************************************************************
 * utility functions
 ************************************************************************/

/* remove trailing path separators */
char *trim_path(char *path)
{
    int last = strlen(path) - 1;
    while ((path[last] == '/') || (path[last] == '\\')) --last;
    path[++last] = '\0';
    return path;
}

/* test if a file exists */
int file_exists(const char *path)
{
#if defined(_WIN32)
    struct _stat s;

    if (path == NULL) return 0;
    if (_stat(path, &s) != 0) return 0;
    return 1;
#else
    struct stat s;

    if (path == NULL) return 0;
    if (stat(path, &s) != 0) return 0;
    return 1;
#endif
}

/* simple integer square root */
static int isqrt(int n)
{
    int b = 0;

    while (n >= 0) {
        n = n - b;
        b = b + 1;
        n = n - b;
    }
    return b - 1;
}

/* determine next available prime number */
static int next_prime(const int num)
{
    int nprime, factor, root;

    /* nprime has to be larger than num and odd */
    nprime = num + (num & 1) + 1;
    /* there is always a prime between n and 2n */
    while (nprime < 2 * num) {
        /* brute force division test on odd factors up to sqrt(nprime) */
        root = isqrt(nprime) + 1;
        for (factor = 3; factor < root; factor += 2) {
            if (nprime % factor == 0) break;
        }
        /* if the loop didn't exit early, we have found a prime */
        if (factor >= root) return nprime;
        nprime += 2;
    }
    return nprime;
}

/* FNV hash function */
static unsigned int hash_func(const void *key)
{
    const unsigned char *p = key;
    unsigned int h = 2166136261;

    if (!p) return 0;

    while (*p) {
        h = (h * 16777619) ^ *p;
        ++p;
    }
    return h;
}

/************************************************************************
 * structs for data structures
 ************************************************************************/

/* linked list */

typedef struct _llnode llnode_t;
struct _llnode {
    const char *key;
    llnode_t *next;
};

typedef struct {
    llnode_t *head;
    llnode_t *tail;
    int count;
} llist_t;

/* set */

typedef struct {
    llnode_t *buckets;
    int nbuckets;
    int count;
} set_t;

/* map */

typedef struct _mapnode mapnode_t;
struct _mapnode {
    const char *key;
    set_t *val;
    mapnode_t *next;
};

typedef struct {
    mapnode_t *buckets;
    int nbuckets;
    int count;
} map_t;

/************************************************************************
 * linked list functions
 ************************************************************************/

/* allocate and initialize linked list */
static llist_t *llist_init()
{
    llist_t *ll;
    ll = (llist_t *)malloc(sizeof(llist_t));
    if (ll != NULL) {
        ll->head = (llnode_t *)malloc(sizeof(llnode_t));
        ll->head->next = NULL;
        ll->head->key = NULL;
        ll->tail = ll->head;
        ll->count = 0;
    }
    return ll;
}

/* destroy linked list and free all associated storage */
static void llist_free(llist_t *ll)
{
    llnode_t *tmp;
    if (ll == NULL) return;

    while (ll->head->next) {
        tmp = ll->head;
        ll->head = ll->head->next;
        free((void *)tmp->key);
        free((void *)tmp);
    }
    free((void *)ll->head);
    free((void *)ll);
}

/* append an item to the end of the linked list */
static void llist_append(llist_t *ll, const char *key)
{
    llnode_t *tmp;
    if ((ll == NULL) || (key == NULL)) return;

    ll->tail->key = my_strdup(key);
    ll->count++;
    tmp = (llnode_t *)malloc(sizeof(llnode_t));
    tmp->key = NULL;
    tmp->next = NULL;
    ll->tail->next = tmp;
    ll->tail = tmp;
}

static int llist_size(llist_t *ll)
{
    if (ll) return ll->count;
    return 0;
}

static void llist_print(llist_t *ll)
{
    llnode_t *tmp;
    if (ll == NULL) return;
    tmp = ll->head;

    if (tmp->next) {
        fputs(tmp->key, stdout);
        tmp = tmp->next;
    }

    while (tmp->next) {
        fputc(':', stdout);
        fputs(tmp->key, stdout);
        tmp = tmp->next;
    }
    fputc('\n', stdout);
}

/************************************************************************
 * set functions
 ************************************************************************/

/* initialize empty set */
static set_t *set_init(int num)
{
    set_t *s = (set_t *)malloc(sizeof(set_t));

    s->nbuckets = next_prime(num);
    s->buckets = malloc(s->nbuckets * sizeof(llnode_t));
    memset(s->buckets, 0, s->nbuckets * sizeof(llnode_t));
    s->count = 0;

    return s;
}

/* destroy a set and free the associated storage */
static void set_free(set_t *s)
{
    llnode_t *tmp, *next;
    int i;

    if (!s) return;
    for (i = 0; i < s->nbuckets; ++i) {
        tmp = s->buckets + i;
        while (tmp->next != NULL) {
            free((void *)tmp->key);
            next = tmp->next->next;
            tmp->key = tmp->next->key;
            free((void *)tmp->next);
            tmp->next = next;
        }
    }
    free((void *)s->buckets);
    free((void *)s);
}

/* add an entry to the set */
static void set_add(set_t *s, const char *key)
{
    llnode_t *tmp;
    unsigned int idx;

    if (!s) return;

    idx = hash_func(key) % s->nbuckets;
    tmp = s->buckets + idx;
    while (tmp->next != NULL) {
        if (strcmp(tmp->key, key) == 0) return;
        tmp = tmp->next;
    }
    s->count++;
    tmp->key = my_strdup(key);
    tmp->next = (llnode_t *)malloc(sizeof(llnode_t));
    tmp = tmp->next;
    tmp->key = NULL;
    tmp->next = NULL;
}

/* find an entry in the set */
static int set_find(set_t *s, const char *key)
{
    llnode_t *tmp;
    unsigned int idx;

    if (!s) return 0;

    idx = hash_func(key) % s->nbuckets;
    tmp = s->buckets + idx;
    while (tmp->next != NULL) {
        if (strcmp(tmp->key, key) == 0) return 1;
        tmp = tmp->next;
    }
    return 0;
}

static int set_size(set_t *s)
{
    if (s) return s->count;
    return 0;
}

/************************************************************************
 * map functions
 ************************************************************************/

/* initialize empty map */
static map_t *map_init(int num)
{
    map_t *m = (map_t *)malloc(sizeof(map_t));
    if (!m) return NULL;

    m->nbuckets = next_prime(num);
    m->buckets = malloc(m->nbuckets * sizeof(mapnode_t));
    memset(m->buckets, 0, m->nbuckets * sizeof(mapnode_t));
    m->count = 0;

    return m;
}

/* destroy a map and free the associated storage */
static void map_free(map_t *m)
{
    mapnode_t *tmp, *next;
    int i;

    if (!m) return;
    for (i = 0; i < m->nbuckets; ++i) {
        tmp = m->buckets + i;
        while (tmp->next != NULL) {
            free((void *)tmp->key);
            set_free(tmp->val);
            next = tmp->next->next;
            tmp->key = tmp->next->key;
            tmp->val = tmp->next->val;
            free((void *)tmp->next);
            tmp->next = next;
        }
    }
    free((void *)m->buckets);
    free((void *)m);
}

/* add an entry to the map */
static void map_add(map_t *m, const char *key, const char *val)
{
    mapnode_t *tmp;
    unsigned int idx;

    if (!m) return;

    idx = hash_func(key) % m->nbuckets;
    tmp = m->buckets + idx;
    while (tmp->next != NULL) {
        if (strcmp(tmp->key, key) == 0) break;
        tmp = tmp->next;
    }

    /* add new entry to map */
    if (tmp->next == NULL) {
        m->count++;
        tmp->key = my_strdup(key);
        tmp->val = set_init(50); /* XXX: chosen arbitrarily */
        tmp->next = (mapnode_t *)malloc(sizeof(mapnode_t));
        tmp->next->key = NULL;
        tmp->next->val = NULL;
        tmp->next->next = NULL;
    }
    set_add(tmp->val, val);
}

/* return an entry in the map */
static set_t *map_find(map_t *m, const char *key)
{
    mapnode_t *tmp;
    unsigned int idx;

    if (!m) return 0;

    idx = hash_func(key) % m->nbuckets;
    tmp = m->buckets + idx;
    while (tmp->next != NULL) {
        if (strcmp(tmp->key, key) == 0) return tmp->val;
        tmp = tmp->next;
    }
    return NULL;
}

static int map_size(map_t *m)
{
    if (m) return m->count;
    return 0;
}

/************************************************************************/

/* combine search for file across paths */
static void make_path(const char *file, llist_t *paths, char *buffer)
{
    llnode_t *tmp;
    const char *val;
    int i;

    tmp = paths->head;
    buffer[0] = '\0';
    while (tmp->next != NULL) {
        i = 0;
        val = tmp->key;

        while (*val) buffer[i++] = *val++;
#if defined(_WIN32)
        buffer[i++] = '\\';
#else
        buffer[i++] = '/';
#endif
        val = file;
        while (*val) buffer[i++] = *val++;
        buffer[i] = '\0';

        if (file_exists(buffer)) return;
        tmp = tmp->next;
    }
    buffer[0] = '\0';
}

/************************************************************************/

static void find_includes(llnode_t *head, llist_t *todo, llist_t *paths,
                          set_t *incl, map_t *deps)
{
    FILE *fp;
    llnode_t *tmp;
    char *buffer, *full, *ptr, *end;
    const char *file;

    buffer = (char *)malloc(4096);
    full = (char *)malloc(4096);

    tmp = head;
    while (tmp->next != NULL) {
        file = tmp->key;
        fp = fopen(file, "r");
        if (fp == NULL) {
            perror("Cannot read source");
            fprintf(stderr, "For file: %s\n", file);
            exit(EXIT_FAILURE);
        }

        /* read file line by line and look for #include "..." */
        while (!feof(fp) && !ferror(fp)) {
            if (fgets(buffer, 4096, fp) == NULL) continue;
            ptr = buffer;
            while (*ptr == ' ' || *ptr == '\t') ++ptr;
            if (*ptr != '#') continue;
            while (*ptr == ' ' || *ptr == '\t') ++ptr;
            if (*++ptr != 'i') continue;
            if (*++ptr != 'n') continue;
            if (*++ptr != 'c') continue;
            if (*++ptr != 'l') continue;
            if (*++ptr != 'u') continue;
            if (*++ptr != 'd') continue;
            if (*++ptr != 'e') continue;
            ++ptr;
            while (*ptr == ' ' || *ptr == '\t') ++ptr;
            if (*ptr != '"') continue;
            ++ptr;
            end = ptr;
            while (*end != '"') {
                if (*end == '\0') {
                    fprintf(stderr, "Unmatched '\"': %s\n", buffer);
                    exit(EXIT_FAILURE);
                }
                ++end;
            }
            *end = '\0';

            /* get full path to include file */
            make_path(ptr, paths, full);
            /* skip, if not found or unreadable. */
            if (full[0] == '\0') continue;

            /* if this is a yet unknown include, add to the
             * todo list, if append is enabled */
            if (set_find(incl, full) == 0) {
                set_add(incl, full);
                llist_append(todo, full);
            }

            map_add(deps, file, full);
        }
        fclose(fp);
        tmp = tmp->next;
    }
    free(buffer);
    free(full);
}

/************************************************************************/

static void add_depend(const char *source, set_t *incl, map_t *deps)
{
    set_t *mydeps;
    llnode_t *tmp;
    int i, num;

    if (source == NULL) return;

    mydeps = map_find(deps, source);
    if (mydeps != NULL) {
        num = mydeps->nbuckets;

        for (i = 0; i < num; ++i) {
            tmp = mydeps->buckets + i;
            while (tmp->next != NULL) {
                if (set_find(incl, tmp->key) == 0) {
                    set_add(incl, tmp->key);
                    add_depend(tmp->key, incl, deps);
                }
                tmp = tmp->next;
            }
        }
    }
}

/************************************************************************/

static void do_depend(llnode_t *head, map_t *deps)
{
    llnode_t *tmp, *lnk;
    set_t *incl;
    const char *source;
    char *target, *ptr;
    int i, num, ext;

    tmp = head;
    while (tmp->next != NULL) {
        source = tmp->key;
        target = strrchr(source, '/');
        if (target == NULL) {
            target = my_strdup(source);
        } else {
            target = my_strdup(target + 1);
        }

        ext = 0;
        ptr = strrchr(target, '.');
        if (ptr != NULL) {
            for (i = 0; i < numextensions; ++i) {
                if (strcmp(ptr, extensions[i]) == 0) ++ext;
            }
            if (ext > 0) {
                ptr[1] = 'o';
                ptr[2] = '\0';
            }
        }

        if (ext > 0) {
            fputs(target, stdout);
            fputs(" : ", stdout);
            fputs(source, stdout);

            incl = set_init(50);
            add_depend(source, incl, deps);

            num = incl->nbuckets;
            for (i = 0; i < num; ++i) {
                lnk = incl->buckets + i;
                while (lnk->next != NULL) {
                    fputc(' ', stdout);
                    fputs(lnk->key, stdout);
                    lnk = lnk->next;
                }
            }
            fputc('\n', stdout);
            set_free(incl);
        }

        free((void *)target);
        tmp = tmp->next;
    }
}

/************************************************************************/

int main(int argc, char **argv)
{
    llist_t *paths, *src, *todo;
    set_t *incl;
    map_t *deps;

    if (argc < 2) {
        fprintf(stderr,
                "FastDep v%s for LAMMPS\n"
                "Usage: %s [-I <path> ...] -- <src1> [<src2> ...]\n",
                version, argv[0]);
        fprintf(stderr, "Supported extensions: %d, %s, %s\n", numextensions,
                extensions[0], extensions[1]);
        return 1;
    }

    /* hash tables for all known included files and dependencies
     * we guesstimate a little over 2x as many entries as sources. */
    incl = set_init(2 * argc);
    deps = map_init(2 * argc);

    /* list of include search paths. prefixed by "." and "..". */
    paths = llist_init();
    llist_append(paths, ".");
    llist_append(paths, "..");

    while (++argv, --argc > 0) {
        if (strncmp(*argv, "-I", 2) == 0) {
            if ((*argv)[2] != '\0') {
                llist_append(paths, trim_path(*argv + 2));
            } else {
                ++argv;
                --argc;
                if (argc > 0) {
                    if (strcmp(*argv, "--") == 0) {
                        break;
                    } else {
                        llist_append(paths, trim_path(*argv));
                    }
                }
            }
        } else if (strcmp(*argv, "--") == 0) {
            break;
        } /* ignore all unrecognized arguments before '--'. */
    }

    src = llist_init();
    while (++argv, --argc > 0) { llist_append(src, *argv); }

    /* process files to look for includes */
    todo = llist_init();
    find_includes(src->head, todo, paths, incl, deps);
    find_includes(todo->head, todo, paths, incl, deps);
    llist_free(todo);

    fprintf(stdout, "# FastDep v%s for LAMMPS\n", version);
    fputs("# Search path: ", stdout);
    llist_print(paths);
    fprintf(stdout, "# % 5d sources\n# % 5d includes\n# % 5d depfiles\n",
            llist_size(src), set_size(incl), map_size(deps));

    set_free(incl);
    do_depend(src->head, deps);

    llist_free(src);
    llist_free(paths);
    map_free(deps);
    return 0;
}

/*
 * Local Variables:
 * compile-command: "gcc -o fastdep.exe -Wall -g -O fastdep.c"
 * c-basic-offset: 4
 * End:
 */
