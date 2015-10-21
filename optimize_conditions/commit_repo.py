import os
import sys
import sh


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print "need to supply message"
        exit()

    rnamake = os.environ['RNAMAKE']
    repo = sh.git.bake(_cwd=rnamake)
    repo.add("-u")
    repo.commit(m=sys.argv[1])
    repo.push("origin")
    print repo.log("-1", '--format="%H"')
