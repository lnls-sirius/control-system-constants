ARG FAC_IMG_HTTPD_TAG=2.4.52-bullseye

FROM httpd:${FAC_IMG_HTTPD_TAG}

WORKDIR /usr/local/apache2/htdocs/

COPY . /usr/local/apache2/htdocs/
