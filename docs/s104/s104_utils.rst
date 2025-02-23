The functions in s104.utils can be accessed by just importing s104.  They are automatically imported
but can be brought in separately if desired.::

    from s100py import s104
    s104.from_gdal(...)
    # or
    from s100py.s104 import utils
    utils.from_gdal(...)
    # or
    from s100py import s104.utils
    s104.utils.from_gdal(...)

..  automodapi:: s100py.s104.v2_0.utils
