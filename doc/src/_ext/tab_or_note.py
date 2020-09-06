
def replace_tabs_handler(app, docname, source):
    """ When builder is not 'html', remove 'tabs' directive
    and replace 'tab' directive with 'admonition'"""
    if app.builder.name != 'html':
        for i in range(len(source)):
            str = source[i].replace('.. tabs::','')
            str = str.replace('.. tab::','.. admonition::')
            source[i] = str

def setup(app):
    app.connect('source-read', replace_tabs_handler)
    return {
        'version': '0.1',
        'parallel_read_safe': True,
        'parallel_write_safe': True,
    }
