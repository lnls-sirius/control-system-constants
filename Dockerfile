ARG FAC_IMG_HTTPD_TAG=2.4.52-bullseye

FROM httpd:${FAC_IMG_HTTPD_TAG}

WORKDIR /usr/local/apache2/htdocs/

COPY ./apache2_conf/httpd.conf /usr/local/apache2/conf/httpd.conf
COPY ./apache2_conf/httpd-mpm.conf /usr/local/apache2/conf/extra/httpd-mpm.conf

COPY . /usr/local/apache2/htdocs/
