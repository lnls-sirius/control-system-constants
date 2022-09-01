#!/usr/bin/env python-sirius

"""
#############################################################################
### Run this script to update 'gop-procedures.html' with the current file ###
### links found in the SharePoint 'Procedimentos' folder. You have to     ###
### have access to the files to run it properly.                          ###
#############################################################################

example:
    build_link_page.py > gop-procedures.html
"""

import sys as _sys
from getpass import getpass as _getpass
from html import entities as _entities
import json as _json
import sharepy as _sharepy


class ClientSharePoint:
    """."""

    CNPEM_URL = 'https://cnpemcamp.sharepoint.com'
    SITE = None  # Example: 'GOP'
    FOLDER = None   # Example: 'Procedimentos'

    _HEADERS = {'accept': 'application/json'}
    _TABLE = None

    def __init__(self):
        """."""
        self._username = None
        self._connector = None
        self._response = None
        self._response_json = None

    def connect(self, username=None, password=None):
        """."""
        if not username:
            stdout = _sys.stdout
            _sys.stdout = _sys.stderr
            print('Username: ', end='')
            _sys.stdout = stdout
            username = input()
        self._username = username
        password = password or _getpass('Password: ')
        # NOTE: manage possible return errors
        self._connector = _sharepy.connect(
            ClientSharePoint.CNPEM_URL, self._username, password)

    def get_response(self):
        """."""
        fullurl = ClientSharePoint.get_fullurl()
        self._response = self._connector.get(
            fullurl, headers=ClientSharePoint._HEADERS)
        self._response_json = ClientSharePoint._convert_to_json(self._response)

    @property
    def files_json(self):
        """."""
        return self._response_json['value']

    @staticmethod
    def get_fullurl():
        """."""
        site = ClientSharePoint.SITE
        folder = 'Documentos%20Compartilhados/' + ClientSharePoint.FOLDER
        rel_path = '/sites/' + site + '/' + folder + '/'
        api = '/_api/web/GetFolderByServerRelativeUrl'
        fullurl = ClientSharePoint.CNPEM_URL + \
            '/sites/' + site + api + \
            "('" + rel_path + "')/Files"
        return fullurl

    @staticmethod
    def _convert_to_json(files):
        text = ''
        for file in files:
            # print('h')
            token = file.decode('utf-8')
            text += (token)
        text_json = _json.loads(text)
        return text_json

    @staticmethod
    def _html_link(file_link, file_name):
        start_link = "&#183; <a href='"
        end_link = "</a>\n<br>\n"
        text_link = start_link + file_link + "'>" + file_name + end_link
        return text_link

    @staticmethod
    def _create_charconv_table():
        items = _entities.codepoint2name.items()
        table_ = {k: '&{};'.format(v) for k, v in items}
        table = dict()
        for key, val in table_.items():
            if 'acute' in val or 'cedil' in val or 'circ' in val or \
              'tilde' in val:
                table[key] = val
        return table


def generate_html(csp):
    """."""
    csp.get_response()

    page_base = (
        "<!DOCTYPE html>\n"
        "<html lang='pt'>\n"
        "<head>\n"
        "</head>\n"
        "<body>\n"
        "<h1>Procedimentos do GOP</h1>\n\r"
        )

    if not ClientSharePoint._TABLE:
        ClientSharePoint._TABLE = ClientSharePoint._create_charconv_table()

    strt = page_base
    for item in csp._response_json['value']:
        file_link = item['LinkingUri']
        file_name = item['Name']
        text_link = ClientSharePoint._html_link(file_link, file_name)
        strt += text_link.translate(ClientSharePoint._TABLE)

    print(strt)


if __name__ == "__main__":

    ClientSharePoint.CNPEM_URL = 'https://cnpemcamp.sharepoint.com'
    ClientSharePoint.SITE = 'GOP'
    ClientSharePoint.FOLDER = 'Procedimentos'

    csp = ClientSharePoint()
    csp.connect(username=None, password=None)
    generate_html(csp)
