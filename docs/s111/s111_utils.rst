The functions in s111.utils can be accessed by just importing s111.  They are automatically imported
but can be brought in separately if desired.::

    from s100py import s111
    s111.from_gdal(...)
    # or
    from s100py.s111 import utils
    utils.from_gdal(...)
    # or
    from s100py import s111.utils
    s111.utils.from_gdal(...)

..  automodapi:: s100py.s111.v1_2.utils
