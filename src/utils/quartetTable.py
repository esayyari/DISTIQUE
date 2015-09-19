import sys
frq=dict()
for line in sys.stdin:
    trip=line.replace("|","").split()
    sortedTrip=sorted(trip)
    k = "/".join(sortedTrip)
    keyt = sorted(line.split(" | "))
    keym= sorted(keyt[0].split(" "))
    
	
    
    v = frq.get(k, dict (zip(trip[1:4],[0]*3)))
    v[keym[1]] += 1
    frq[k] = v

for k in sorted(frq.keys()):
    print k,
    tf = frq[k]
    for t in sorted(tf.keys()):
        print tf[t],
    print
