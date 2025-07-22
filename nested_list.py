# the following class allows easy access of nested lists
# falcon trees are binary trees which are stored as nested lists.

class nestedList(list):
    def __init__(self, *args):
        if len(args) == 1:
            super().__init__(nestedList(arg) if type(arg) is list else arg for arg in args[0])
        else:
            super().__init__(nestedList(arg) if type(arg) is list else arg for arg in args)

    def __getitem__(self, key):
        if type(key) is int:
            return super().__getitem__(key)
        if len(key) == 1:
            return super().__getitem__(key[0])
        return self.__getitem__(key[0])[key[1:]]

    def __setitem__(self, key, value):
        if type(key) is int:
            super().__setitem__(key, value)
            return
        if len(key) == 1:
            super().__setitem__(key[0], value)
            return
        self.__getitem__(key[0]).__setitem__(key[1:], value)
