#############################################################################
### Run this script to update 'gop-procedures.html' with the current file ###
### links found in the SharePoint 'Procedimentos' folder. You have to     ###
### have access to the files to run it properly.                          ###
#############################################################################

from pwinput import pwinput
import sharepy, os, json

site_url = 'https://cnpemcamp.sharepoint.com'
username = input('Username: ')
passw = pwinput('Enter password: ')

s = sharepy.connect(site_url, username, passw)

headers = {'accept': 'application/json'}
rel_path = '/sites/GOP/Documentos%20Compartilhados/Procedimentos/'
api = '/_api/web/GetFolderByServerRelativeUrl'
fullurl = site_url + '/sites/GOP' + api + "('" + rel_path + "')/Files"
print(fullurl)

files = s.get(fullurl, headers=headers)

dir_path = os.path.dirname(os.path.realpath(__file__))
file_path = os.path.join(dir_path, 'gop-procedures.html')
r  = ''

page_base = "<!DOCTYPE html>\n\
<html lang='pt'>\n\
<head>\n\
</head>\n\
<body>\n\
<h1>Procedimentos do GOP</h1>\n\r"

start_link = "&#183; <a href='"
end_link = "</a>\n<br>\n"

for file in files:
    r += (file.decode('utf-8'))
r_json = json.loads(r)

open(file_path, 'w').close()
f = open(file_path, 'a')
f.write(page_base)
for item in r_json['value']:
    file_link = item['LinkingUri']
    file_name = item['Name']
    f.write(start_link + file_link + "'>" + file_name + end_link)
f.close()

