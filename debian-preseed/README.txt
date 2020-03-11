Debian Preseed Files
====================

Reference: https://wiki.debian.org/DebianInstaller/Preseed

To quote the Debian page:

    Preseeding provides a way to set answers to questions
    asked during the installation process, without having
    to manually enter the answers while the installation
    is running. This makes it possible to fully automate
    most types of installation and even offers some
    features not available during normal installations.

There are numerous ways to provide a preseed file to the
installation process. One of the most simple ones is to supply
a URL to the boot parameters prior to the installer.

In this case, after the graphical Debian installer opens up,

  1. Select the "Graphical Install",
  2. then "Graphical Automated Install..." option.
  3. Inside this menu, highlight the "Automated Install" option but do NOT press ENTER.
     Instead, press the letter "e" and change the boot parameters.

     The boot parameters will be something like:

         auto=true priority=critical vga=788 initrd=/install.amd/initrd.gz --- quiet

     So, add the following option prior to "---":

         url=http://<path/to/preseed/file>

     The final line should be something like:

         auto=true priority=critical vga=788 initrd=/install.amd/initrd.gz url=http://<path/to/preseed/file> ---  quiet

  4. After editing the line press "F10" to start the installation process.
