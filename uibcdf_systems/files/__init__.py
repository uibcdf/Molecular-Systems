from os import listdir, path

cdir = path.dirname(path.abspath(__file__))

for dirname in sorted(listdir(cdir)):
    absdirname = path.join(cdir, dirname)
    if path.isdir(absdirname):
        dict_files = {}
        for filename in sorted(listdir(absdirname)):
            absfilename = path.join(absdirname, filename)
            if path.isfile(absfilename) and not absfilename.endswith('.py'):
                dict_files[filename]=absfilename
        globals()[dirname]=dict_files

del(cdir, dict_files, listdir, path, dirname, absdirname, filename, absfilename)
