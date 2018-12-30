class DataFile:
    def open(cls, path, mode):
        """
        Open a data file.  The base implementation of this function just
        opens and returns a standard python file object.
        Subclasses can override to either open files using different openers
        (like fitsio.FITS), or, for more specific data types, return an
        instance of the class itself to use as an intermediary for the file.
        """
        return open(path, mode)

class FitsFile(DataFile):
    """
    A data file in the FITS format.
    Using these files requires the fitsio package.
    """
    suffix = 'fits'
    
    @classmethod
    def open(cls, path, mode, **kwargs):
        import fitsio
        # Fitsio doesn't have pure 'w' modes, just 'rw'.
        # Maybe we should check if the file already exists here?
        if mode == 'w':
            mode = 'rw'
        return fitsio.FITS(path, mode=mode, **kwargs)
