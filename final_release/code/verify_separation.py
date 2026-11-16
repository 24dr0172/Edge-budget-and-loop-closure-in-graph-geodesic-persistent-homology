"""Cyclic orientation search for finite threshold graphs."""
def cyclic_orientation(order,edges):
    """Solve the orientation implications; return a full directed certificate.

    Each possible direction is a literal. It implies both required directions
    through every intervening vertex. A missing edge forces its negation.
    Transitive closures plus exhaustive branching over unassigned edges give
    satisfiability; every returned assignment is checked directly.
    """
    literals=[d for a,b in edges for d in [(a,b),(b,a)]]
    ix={d:i for i,d in enumerate(literals)}
    reach=[1<<i for i in range(len(literals))]
    pos={v:i for i,v in enumerate(order)}
    for i,(a,b) in enumerate(literals):
        for step in range(1,(pos[b]-pos[a])%6):
            w=order[(pos[a]+step)%6]
            for pair in [(a,w),(w,b)]:
                reach[i]|=1<<(ix[pair] if pair in ix else i^1)
    for k in range(len(literals)):
        for i in range(len(literals)):
            if reach[i]>>k&1:reach[i]|=reach[k]
    even=sum(1<<i for i in range(0,len(literals),2));odd=even<<1
    def solve(mask):
        if mask&(((mask&even)<<1)|((mask&odd)>>1)):return None
        for i in range(0,len(literals),2):
            if not mask&(3<<i):
                a=solve(mask|reach[i])
                return a if a is not None else solve(mask|reach[i+1])
        return mask
    mask=solve(0)
    if mask is None:return None
    directed=[pair for i,pair in enumerate(literals) if mask>>i&1]
    assert orientation_valid(order,set(directed))
    return directed


def orientation_valid(order,directed):
    pos={v:i for i,v in enumerate(order)}
    return all(all((a,order[(pos[a]+step)%6]) in directed and
                   (order[(pos[a]+step)%6],b) in directed
                   for step in range(1,(pos[b]-pos[a])%6))
               for a,b in directed)


