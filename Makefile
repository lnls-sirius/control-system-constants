install-html:
	mkdir -p /var/www/html
	cp -ra * /var/www/html/
	find /var/www/html -type d -exec chmod +x {} \;
