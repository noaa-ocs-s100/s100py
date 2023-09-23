The classes in s111.api can be accessed by just importing s111.  They are automatically imported
but can be brought in separately if desired.::

    from s100py import s111
    s111.s111File(...)
    # or
    from s100py.s111 import api
    api.s111File(...)
    # or
    from s100py import s111.api
    s111.api.s111File(...)

..  automodapi:: s100py.s111.v1_2.api
