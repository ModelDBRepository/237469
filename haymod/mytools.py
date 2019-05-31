import math

def spike_times(time,vrec,V_min_peak=-20,V_max_valley=0):
  valley_reached = 1
  sptime = []  
  for j in range(1,len(time)-1):
    if valley_reached and vrec[j] >= V_min_peak and vrec[j] > vrec[j-1] and vrec[j] >= vrec[j+1]:
      valley_reached = 0
      sptime.append(time[j])
    elif valley_reached==False and vrec[j] <= V_max_valley:
      valley_reached = 1
  return sptime  

#membpotderivs(time,vrec): Given the membrane potentials (vrec) at time points time[0],time[1],...,time[N],
#return the derivatives at time points time[1],time[2],...,time[N-1]
def membpotderivs(time,vrec):
  N = len(time)
  tdiff = [x-y for x,y in zip(time[1:N-1],time[0:N-2])]
  vdiff = [x-y for x,y in zip(vrec[1:N-1],vrec[0:N-2])]
  mderiv = [x/y for x,y in zip(vdiff,tdiff)]
  return [0.5*(x+y) for x,y in zip(mderiv[1:N-2],mderiv[0:N-3])]

#limitcyclescaledv(v1,dv1,v2,dv2): Give the coefficient for memb. pot. derivative that one has to use in order to make
#the difference on the derivative axis as significant as the difference on the memb. pot. axis
def limitcyclescaledv(v1,dv1,v2,dv2):
  maxv = max(max(v1),max(v2))
  minv = min(min(v1),min(v2))
  maxdv = max(max(dv1),max(dv2))
  mindv = min(min(dv1),min(dv2))
  return 1.0*(maxv-minv)/(maxdv-mindv)

def limitcyclediff(v1,dv1,v2,dv2,dvcoeff=0.1):
  N1 = len(v1)
  N2 = len(v2)
  dv1 = [dvcoeff*x for x in dv1]
  dv2 = [dvcoeff*x for x in dv2]
  Dmin = N1*[0]
  for i in range(0,N1):
    Dmin[i] = math.sqrt(min([(x-v1[i])**2+(y-dv1[i])**2 for x,y in zip(v2,dv2)]))
  vdiff = [x-y for x,y in zip(v1[1:N1],v1[0:N1-1])] 
  dvdiff = [x-y for x,y in zip(dv1[1:N1],dv1[0:N1-1])] 
  h = [math.sqrt(x**2+y**2) for x,y in zip(vdiff,dvdiff)]
  #use the trapezoid rule for integration:
  Dminmean = [(x+y)/2.0 for x,y in zip(Dmin[1:N1],Dmin[0:N1-1])]
  print "hsum="+str(sum(h))
  return sum([x*y for x,y in zip(Dminmean,h)])
  

def interpolate(tref,vref,tint): #Assumes that the trefs come sorted!
  vint = len(tint)*[0.0]
  addedOne = False
  #print tref
  #print tint
  #if tref[len(tref)-1] == tint[len(tint)-1]:
  #  tref.append(tref[len(tref)-1]+0.0001)
  #  vref.append(vref[len(tref)-1])
  #  addedOne = True
  if tref[0] > tint[0] or tref[len(tref)-1] < tint[len(tint)-1]:
    print "Extrapolation needed!"
    return len(tint)*[-1]
  indvrecnow = 0  
  for j in range(0,len(tint)):
    while tref[indvrecnow+1] <= tint[j]:
      indvrecnow = indvrecnow + 1
      if indvrecnow == len(tref)-1: # It must be the last index if this happens
        vint[j:len(tint)] = [vref[indvrecnow]]*(len(tint)-j)
        return vint
    vint[j] = vref[indvrecnow] + 1.0*(tint[j]-tref[indvrecnow])/(tref[indvrecnow+1]-tref[indvrecnow])*(vref[indvrecnow+1]-vref[indvrecnow])
  return vint

#kronecker product of list A and list B
def kron(A,B):
  C = []
  if type(B[0]) is int or type(B[0]) is float:
    for i in range(0,len(A)):
      for j in range(0,len(B)):
        print "asdf"
        print B[j]
        C.append(A[i]*B[j])
  elif type(B[0][0]) is int or type(B[0][0]) is float:
    for i in range(0,len(A)):
      for j in range(0,len(B)):
        C.append([x*A[i] for x in B[j]])
  return C
  
def cumprod(A):
  B = len(A)*[0]; B[0]=A[0]
  for j in range(1,len(A)):
    B[j] = B[j-1]*A[j]
  return B
