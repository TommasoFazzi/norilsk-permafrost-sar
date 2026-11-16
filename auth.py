"""Earthdata Login credentials, read from ~/.netrc or environment.

Never hard-code credentials and never commit .netrc / .env to git.
~/.netrc line:  machine urs.earthdata.nasa.gov login USER password PASS
(chmod 600 ~/.netrc)
"""
import netrc
import os
from getpass import getpass

HOST = "urs.earthdata.nasa.gov"


def earthdata_credentials() -> tuple[str, str]:
    user = os.environ.get("EARTHDATA_USERNAME")
    pwd = os.environ.get("EARTHDATA_PASSWORD")
    if user and pwd:
        return user, pwd
    try:
        auth = netrc.netrc().authenticators(HOST)
        if auth:
            return auth[0], auth[2]
    except (FileNotFoundError, netrc.NetrcParseError):
        pass
    return input("Earthdata username: "), getpass("Earthdata password: ")


def asf_session():
    import asf_search as asf
    user, pwd = earthdata_credentials()
    return asf.ASFSession().auth_with_creds(user, pwd)
