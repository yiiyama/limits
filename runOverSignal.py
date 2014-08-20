import sys
import os
import subprocess
import time

thisdir = os.path.dirname(os.path.abspath(__file__))

if __name__ == '__main__':

    model = sys.argv[1]
    script = sys.argv[2]
    nodePool = sys.argv[3]

    if model == 'TChiwg':
        points = ['TChiwg_%d' % mchi for mchi in range(200, 810, 10)]
    elif model == 'T5wg':
        points = ['T5wg_%d_%d' % (mglu, mchi) for mglu in range(700, 1350, 50) for mchi in range(25, mglu, 50)]

    procs = []
    while True:
        running = filter(lambda x: x[1].poll() is None, procs)
        done = list(set(procs) - set(running))
        
        for point, proc, node in done:
            if proc.poll() == 0:
                print point, 'Done'
            else:
                print '!!!!!!!', point, 'Return Code', proc.poll()

        procs = running

        if len(procs) >= 16:
            time.sleep(10)

        if len(points) == 0:
            if len(procs) != 0:
                print [x[0] for x in procs]
                time.sleep(10)
                continue
            else:
                break

        usedNodes = [x[2] for x in procs]

        hostProc = subprocess.Popen(['host', nodePool], stdout = subprocess.PIPE, stderr = subprocess.PIPE)
        out, err = hostProc.communicate()
        for line in out.split('\n'):
            if 'has address' not in line: continue
            addr = line.split()[3]
            if addr not in usedNodes:
                node = addr
                break
        else:
            time.sleep(2)
            continue

        point = points.pop()
        print 'Submit', point, 'to', node
        
        proc = subprocess.Popen('ssh -oStrictHostKeyChecking=no -oLogLevel=quiet %s "screen -D -m %s %s"' % (node, script, point), shell = True)
        procs.append((point, proc, node))
