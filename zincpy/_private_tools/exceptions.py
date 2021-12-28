# Exceptions derived from value errors

class ZincPyValueError(ValueError):
    pass

class InvalidZincIdError(ZincPyValueError):
    pass

class InvalidCatalogError(ZincPyValueError):
    pass

class InvalidFileFormatError(ZincPyValueError):
    pass

class InvalidAvailabilityError(ZincPyValueError):
    pass

class InvalidBioactiveError(ZincPyValueError):
    pass

class InvalidBiogenicError(ZincPyValueError):
    pass

class InvalidReactivityError(ZincPyValueError):
    pass

class NegativeCountError(ZincPyValueError):
    pass

class InvalidSubsetError(ZincPyValueError):
    pass

class InvalidMolecularWeightRangeError(ZincPyValueError):
    pass

class InvalidLogPRangeError(ZincPyValueError):
    pass

class InvaludUrlTypeError(ZincPyValueError):
    pass

# Exceptions derived from type errors

class ZincPyTypeError(TypeError):
    pass

class CountTypeError(ZincPyTypeError):
    pass

# IO Errors

class ZincPyIOError(IOError):
    pass

class DownloadError(ZincPyIOError):
    pass

class ZincNotFoundError(ZincPyIOError):
    pass

class ZincTimeoutError(ZincPyIOError):
    pass