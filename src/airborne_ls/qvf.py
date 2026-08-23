"""
Mimic the functionality of the JRSRP qvf module, for manipulating QVF-style
file names.

"""
import os


def getfields(fullpath):
    """
    Return the fields list of the given full path (reverse of assemblefields)

    Removes directory spec and parts after '.', splits on '_'

    Returns:
      fields (list[str]): Separated fields of the filename
    """
    fn = os.path.basename(fullpath)
    basefn = fn.split('.')[0]
    fields = basefn.split('_')
    if len(fields) < 4:
        msg = f"Filename {fullpath} has only {len(fields)} fields"
        raise ValueError(msg)
    return fields


def assemblefields(fields):
    """
    Assemble the given fields to a base name (reverse of getfields)
    """
    if len(fields) < 4:
        msg = f"Fields list {fields} has only {len(fields)} fields"
        raise ValueError(msg)
    return '_'.join(fields)


def _changebasefn(fullpath, newbase):
    """
    Change the base part of the given full path
    """
    (dirname, filename) = os.path.split(fullpath)
    dotndx = filename.find('.')
    if dotndx != -1:
        suffix = filename[dotndx:]
    else:
        suffix = ""
    return os.path.join(dirname, newbase + suffix)


def getstagecode(filename):
    """
    Get the processing stage code from the given file name
    """
    fields = getfields(filename)
    stagecode = fields[3][:3]
    return stagecode


def setstagecode(filename, stage):
    """
    Change the stage code of the given filename

    Parameters:
      filename (str): Filename to change
      stage (str): New stage code

    Returns:
      newfile (str): Filename with new stage code
    """
    fields = getfields(filename)
    stageAndZone = fields[3]
    zone = stageAndZone[3:]
    newStageAndZone = stage + zone
    fields[3] = newStageAndZone
    newbase = assemblefields(fields)
    newFilename = _changebasefn(filename, newbase)
    return newFilename


def getoptionfield(fullpath, tagChar):
    """
    Return the value of the tagChar' option field (None if not present)
    """
    fields = getfields(fullpath)
    optDict = {f[0]: f[1:] for f in fields[4:]}
    val = optDict.get(tagChar)
    return val


def setoptionfield(fullpath, tagChar, fieldVal):
    """
    Change the value of the option field for the given tag. If this tag
    is not already present, it will be added. If the fieldVal is None,
    the tag field will be removed.

    Return the new fullpath
    """
    fields = getfields(fullpath)
    optDict = {f[0]: f[1:] for f in fields[4:]}
    if fieldVal is not None:
        optDict[tagChar] = fieldVal
    elif tagChar in optDict:
        optDict.pop(tagChar)
    tagList = sorted(optDict.keys())
    newFields = fields[:4] + [(t + optDict[t]) for t in tagList]
    newbase = assemblefields(newFields)
    newFullpath = _changebasefn(fullpath, newbase)
    return newFullpath


def getwhere(fullpath):
    """
    Return the 'where' field of the given filename
    """
    where = getfields(fullpath)[1]
    return where


def setwhere(fullpath, where):
    """
    Change the 'where' field of the given fullpath to the given value. Return a
    new fullpath.
    """
    fields = getfields(fullpath)
    fields[1] = where
    newbase = assemblefields(fields)
    newFullpath = _changebasefn(fullpath, newbase)
    return newFullpath


def getutmzone(filename):
    """
    Return the UTM zone from the zoneCode field. Returns None
    if not UTM
    """
    fields = getfields(filename)
    zoneField = fields[3][3:]
    if zoneField[0] == 'm':
        utmZone = int(zoneField[1:])
        if len(zoneField) == 2:
            utmZone = utmZone + 50
    else:
        utmZone = None
    return utmZone


def setsuffix(fullpath, suffix):
    """
    Change the file suffix (anything after the first dot in the file name) to
    the given new suffix
    """
    (dirname, filename) = os.path.split(fullpath)
    dotndx = filename.find('.')
    if len(suffix) > 0:
        suffixWithDot = f".{suffix}"
    else:
        suffixWithDot = ''
    if dotndx != -1:
        newFilename = filename[:dotndx] + suffixWithDot
    else:
        newFilename = filename + suffixWithDot
    return os.path.join(dirname, newFilename)


def makeTileWhere(x, y, utmZone):
    """
    Make the 'where' field for a tile. Coordinates are in UTM.

    Parameters:
      x, f (float): The X & Y coordinates of the top-left corner of the tile
      zone (int): UTM zone number of the coordinates

    Returns:
      where (str): The 'where' field for the tile
    """
    xstr = "{:06}".format(x)
    northSouthCode = "n"
    if utmZone < 0:
        northSouthCode = "s"
    ystr = "{}{:07d}".format(northSouthCode, y)
    where = f"x{xstr}y{ystr}z{abs(utmZone)}"
    return where


def makeProjectionCode(utmZone):
    """
    Make a QVF-style projection code for the given UTM zone
    """
    if utmZone != 0:
        zoneNum = abs(utmZone)
        if 50 <= zoneNum <= 56:
            zoneNum = zoneNum - 50
        projCode = f"m{zoneNum}"
    return projCode
