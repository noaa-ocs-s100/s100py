The classes in s104.api can be accessed by just importing s104.  They are automatically imported
but can be brought in separately if desired.::

    from s100py import s104
    s104.s104File(...)
    # or
    from s100py.s104 import api
    api.s104File(...)
    # or
    from s100py import s104.api
    s104.api.s104File(...)

..  automodapi:: s100py.s104.v2_0.api
