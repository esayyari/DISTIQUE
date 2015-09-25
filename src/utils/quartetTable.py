import sys
import re
frq=dict()
for line in sys.stdin:
	posPat=re.match("alpha|gamma|beta",line) 
    	if posPat is not None:
	    line=re.sub("alpha |gamma |beta ","",line)
	    trip=line.replace("|","").split()
	    sortedTrip=sorted(trip)
	    k = "/".join(sortedTrip)
	    keyt = sorted(line.split(" | "))
	    keym= sorted(keyt[0].split(" "))
	    
		
	    
	    v = frq.get(k, dict (zip(trip[1:4],[0.5]*3)))
	    v[keym[1]] += 1
	    frq[k] = v
	else:
	    line = re.sub("delta |epsilon ","",line)
            k = "/".join(sorted(line.split()))
	    trip = sorted(k.split("/"))
            v = frq.get(k, dict (zip(trip[1:4],[0.5]*3)))
            v[trip[1]] += 1
	    v[trip[2]] += 1
	    v[trip[3]] += 1
            frq[k] = v


for k in sorted(frq.keys()):
    print k,
    tf = frq[k]
    for t in sorted(tf.keys()):
	print tf[t],
    print

