sh createdocu.sh
ssh -t mihrke,libeegtools@shell.sourceforge.net create
ssh mihrke@shell.sf.net rm -rf /home/groups/l/li/libeegtools/htdocs/api/current
scp -r doc/html mihrke@web.sf.net:/home/groups/l/li/libeegtools/htdocs/api/current
ssh mihrke@shell.sf.net chmod +r /home/groups/l/li/libeegtools/htdocs/api/current/*


