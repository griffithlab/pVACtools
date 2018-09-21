def compare(path1, path2):
    r1 = open(path1)
    r2 = open(path2)
    result = not len(set(r1.readlines())^set(r2.readlines()))
    r1.close()
    r2.close()
    return result
