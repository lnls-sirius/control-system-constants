ARG FAC_IMG_HTTPD_TAG=2.4.52-bullseye

FROM httpd:${FAC_IMG_HTTPD_TAG}

SHELL ["bash", "-c"]

ARG RELEASE_TAG
RUN echo "RELEASE_TAG=${RELEASE_TAG}"

WORKDIR /usr/local/apache2/htdocs/

RUN set -e; \
    set -x; \
    \
    apt update; \
    apt -y install git nano

COPY ./apache2_conf/httpd.conf /usr/local/apache2/conf/httpd.conf
COPY ./apache2_conf/httpd-mpm.conf /usr/local/apache2/conf/extra/httpd-mpm.conf

#RUN mkdir /usr/local/apache2/htdocs/control-system-constants/
#COPY . /usr/local/apache2/htdocs/control-system-constants/

RUN git clone https://github.com/lnls-sirius/control-system-constants; \
    cd control-system-constants; \
    git checkout ${RELEASE_TAG}
