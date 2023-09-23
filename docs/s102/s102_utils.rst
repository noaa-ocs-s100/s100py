The functions in s102.utils can be accessed by just importing s102.  They are automatically imported
but can be brought in separately if desired.::

    from s100py import s102
    s102.from_gdal(...)
    # or
    from s100py.s102 import utils
    utils.from_gdal(...)
    # or
    from s100py import s102.utils
    s102.utils.from_gdal(...)

..  automodapi:: s100py.s102.v2_2.utils
