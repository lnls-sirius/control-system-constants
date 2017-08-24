git-pull:
	git pull; sudo chown -R fac.fac *
	sudo find ./ -type d -exec chmod +x {} \;

