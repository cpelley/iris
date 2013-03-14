uxterm -bg white -e 'cd ~/git/iris/docs/iris/; make doctest; /bin/bash' &
uxterm -bg white -e 'cd ~/git/iris/docs/iris/; make extest; /bin/bash' &
uxterm -bg white -e 'cd ~/git/iris/; python setup.py test; /bin/bash' &
