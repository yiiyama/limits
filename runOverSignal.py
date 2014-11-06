import sys
import os
import subprocess
import time

def runOverSignal(points, script, nodePool, additionalArgs = []):

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

        model, point = points.pop()
        print 'Submit', model, point, 'to', node
        
        proc = subprocess.Popen('ssh -oStrictHostKeyChecking=no -oLogLevel=quiet {node} "screen -D -m {script} {model} {point} {add}"'.format(node = node, script = script, model = model, point = point, add = ' '.join(additionalArgs)), shell = True)
        procs.append((point, proc, node))


if __name__ == '__main__':

    sPoints = sys.argv[1]
    script = sys.argv[2]
    nodePool = sys.argv[3]

    if ',' in sPoints:
        # model1,point1 model2,point2 ...
        pairs = sPoints.split()
        points = []
        for pair in pairs:
            points.append((pair.partition(',')[0], pair.partition(',')[1]))

    else:
        points = []
        if sPoints == 'All' or sPoints == 'TChiwg':
            points += [('TChiwg', '%d' % mchi) for mchi in range(200, 810, 10)]
        if sPoints == 'All' or sPoints == 'T5wg':
            points += [('T5wg', '%d_%d' % (mglu, mchi)) for mglu in range(700, 1550, 50) for mchi in range(25, mglu, 50)]
        if sPoints == 'All' or sPoints == 'T5wg+TChiwg':
            points += [('T5wg+TChiwg', '%d_%d' % (mglu, mchi)) for mglu in range(700, 1550, 50) for mchi in range(225, mglu, 50)]
        if sPoints == 'All' or sPoints == 'Spectra_gW':
            points += [('Spectra_gW', 'M3_%d_M2_%d' % (m3, m2)) for m3 in range(715, 1565, 50) for m2 in range(205, m3, 50)]

    runOverSignal(points, script, nodePool)
