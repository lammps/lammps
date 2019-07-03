class PathFilter(Filter):
    r"""Filter skipping over simple file paths
    """
    _DOC_ERRORS = ["zA"]
    _pattern = re.compile(r"\/[^\s].*")
    def _skip(self,word):
        if self._pattern.match(word):
            return True
        return False
