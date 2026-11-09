"""Check the finite intrinsic-circle formula by exact boundary reduction.

No persistence package or floating-point distances are used. For each rational
sample, check every pair-distance threshold over F2, F3 and F5, including the
rank of the triangle-boundary span and membership of the consecutive cycle.
These finite checks corroborate the statement; they are not a general proof.
"""
from fractions import Fraction as Q
from itertools import combinations
import json
import random


def reduce_column(column, basis, prime):
    column={i:v % prime for i,v in column.items() if v % prime}
    while column:
        pivot=min(column)
        if pivot not in basis:
            break
        scale=column[pivot]
        for i,value in basis[pivot].items():
            new=(column.get(i,0)-scale*value) % prime
            if new:column[i]=new
            else:column.pop(i,None)
    return column


def add_column(column, basis, prime):
    reduced=reduce_column(column,basis,prime)
    if reduced:
        pivot=min(reduced)
        inverse=pow(reduced[pivot],-1,prime)
        basis[pivot]={i:value*inverse % prime for i,value in reduced.items()}


def check_sample(points):
    points=sorted(points);n=len(points)
    gaps=[points[i+1]-points[i] for i in range(n-1)]+[1+points[0]-points[-1]]
    g=max(gaps)
    assert 0<g<Q(1,2)
    pairs=list(combinations(range(n),2));edge_index={edge:i for i,edge in enumerate(pairs)}
    distance={edge:min(points[edge[1]]-points[edge[0]],1-points[edge[1]]+points[edge[0]]) for edge in pairs}
    triples=list(combinations(range(n),3));winding=[]
    for i,j,k in triples:
        arcs=(points[j]-points[i],points[k]-points[j],1+points[i]-points[k])
        if max(arcs)<Q(1,2):winding.append(max(arcs))
    T=min([Q(1,2)]+winding)
    assert max(g,Q(1,3))<=T<=Q(1,2)
    if all(x==Q(1,n) for x in gaps) and n>=4:
        assert T==Q((n+2)//3,n)
    triangles={}
    for i,j,k in triples:
        threshold=max(distance[i,j],distance[i,k],distance[j,k])
        triangles.setdefault(threshold,[]).append({edge_index[j,k]:1,edge_index[i,k]:-1,edge_index[i,j]:1})
    ring={edge_index[i,i+1]:1 for i in range(n-1)};ring[edge_index[0,n-1]]=-1
    tests=0
    for prime in (2,3,5):
        parent=list(range(n));components=n;basis={};edges=0
        def find(i):
            while parent[i]!=i:
                parent[i]=parent[parent[i]];i=parent[i]
            return i
        for threshold in sorted(set(distance.values()) | {Q(0),Q(1,2)}):
            for (i,j),value in distance.items():
                if value==threshold:
                    edges+=1;a,b=find(i),find(j)
                    if a!=b:parent[a]=b;components-=1
            for column in triangles.get(threshold,[]):add_column(column,basis,prime)
            beta=edges-n+components-len(basis)
            assert beta==int(g<=threshold<T),(points,prime,threshold,g,T,beta)
            if threshold>=g:
                assert bool(reduce_column(ring,basis,prime))==(threshold<T)
            tests+=1
    return tests,g>=Q(1,3),T==g


def verify():
    rng=random.Random(20261028)
    samples=[]
    for n in range(3,11):
        for _ in range(20):
            while True:
                points=[Q(i,360) for i in sorted(rng.sample(range(360),n))]
                gaps=[points[i+1]-points[i] for i in range(n-1)]+[1+points[0]-points[-1]]
                if max(gaps)<Q(1,2):break
            samples.append(points)
    samples.extend([[Q(i,n) for i in range(n)] for n in range(3,31)])
    samples.extend([[Q(0),Q(49,100),Q(1,2),Q(99,100)],
                    list(map(Q,['0','0.08','0.15','0.161','0.29','0.60','0.66']))])
    checks=large=empty=0
    for sample in samples:
        count,is_large,is_empty=check_sample(sample)
        checks+=count;large+=is_large;empty+=is_empty
    result=dict(status='passed',samples=len(samples),threshold_field_checks=checks,
                large_gap_samples=large,empty_reference_samples=empty,coefficient_primes=[2,3,5])
    return result


if __name__=='__main__':print(json.dumps(verify(),indent=2))
