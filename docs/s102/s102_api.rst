The classes in s102.api can be accessed by just importing s102.  They are automatically imported
but can be brought in separately if desired.::

    from s100py import s102
    s102.S102File(...)
    # or
    from s100py.s102 import api
    api.S102File(...)
    # or
    from s100py import s102.api
    s102.api.S102File(...)

..  automodapi:: s100py.s102.v2_2.api
