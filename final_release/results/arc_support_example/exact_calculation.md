# An arc-supported graph whose Rips graph is not cyclic in the curve order

**Yes: six points on a round circle suffice.** Take the consecutive ring, add two chords from one vertex, and use Euclidean chord lengths as edge costs. Every graph edge is arc-supported. Nevertheless, at Rips threshold \(t=19/10\), the threshold graph admits no cyclic orientation in the geometric curve order. Its loop is still present at that threshold.

## 1. The two conditions

Write \(c(x,y)=\|x-y\|\), and let \(D_G\) be the shortest-path metric of a graph \(G\) containing the consecutive ring.

An edge \(\{a,b\}\in E(G)\) is **arc-supported** if one of its curve arcs, listed as

\[
a=w_0,w_1,\ldots,w_m=b,
\]

has an endpoint \(p\in\{a,b\}\) satisfying

\[
c(w_{i-1},w_i)\le c(a,b),\qquad
D_G(p,w_i)\le c(a,b).
\]

The second condition asks for distances from one chosen endpoint to the arc vertices.

For **cyclicity**, an orientation must satisfy

\[
x\prec y\prec z\prec x,\quad x\longrightarrow z
\quad\Longrightarrow\quad
x\longrightarrow y\ \text{and}\ y\longrightarrow z.
\tag{C}
\]

Each undirected edge receives one direction. This is the cyclic-graph convention in [Adams, Coldren and Willmot, Section 3](https://arxiv.org/pdf/1812.03374). Sections 1–7 use the geometric order of the sampled points; Section 8 considers arbitrary cyclic orders throughout the filtration.

## 2. Coordinates and graph

Work first on the circle of radius one. Choose

\[
\begin{aligned}
a&=(1,0), & u&=(24/25,7/25), & v&=(0,1),\\
b&=(-4/5,3/5), & q&=(-1,0), & z&=(0,-1).
\end{aligned}
\]

All six points lie on the circle, in counterclockwise order

\[
a\prec u\prec v\prec b\prec q\prec z\prec a.
\]

Let

\[
E(G)=\{au,uv,vb,bq,qz,za\}\cup\{av,ab\}.
\]

Every edge costs its Euclidean chord length. The relevant lengths are

\[
\begin{aligned}
c(a,u)&=\frac{\sqrt2}{5}, & c(u,v)&=\frac65,
& c(v,b)&=\frac{2}{\sqrt5},\\
c(b,q)&=\frac{\sqrt{10}}5, & c(q,z)&=c(z,a)=\sqrt2,
& c(a,v)&=\sqrt2,\\
c(a,b)&=\frac{3\sqrt{10}}5.
\end{aligned}
\]

For every edge present in \(G\), its direct chord realizes the graph distance: the Euclidean triangle inequality makes every competing path at least as expensive.

## 3. Every edge is arc-supported

Each consecutive ring edge uses its one-edge arc, so it passes immediately.

For the extra edge \(av\), use the arc \(a,u,v\) and anchor \(p=a\). Its consecutive costs are

\[
\frac{\sqrt2}{5},\ \frac65\ \le\ \sqrt2=c(a,v),
\]

and its anchor distances are

\[
D_G(a,a)=0,\quad D_G(a,u)=\frac{\sqrt2}{5},\quad
D_G(a,v)=\sqrt2.
\]

For the extra edge \(ab\), use \(a,u,v,b\), again anchored at \(a\). Its consecutive costs are bounded by

\[
\max\left\{\frac{\sqrt2}{5},\frac65,\frac{2}{\sqrt5}\right\}
=\frac65<\frac{3\sqrt{10}}5=c(a,b).
\]

The anchor distances are

\[
0,\quad\frac{\sqrt2}{5},\quad\sqrt2,\quad\frac{3\sqrt{10}}5,
\]

all at most \(c(a,b)\). Thus **all eight edges of \(G\) are arc-supported**.

## 4. A Rips edge that cannot be oriented either way

Let \(H_t\) be the undirected Rips threshold graph:

\[
\{x,y\}\in E(H_t)\quad\Longleftrightarrow\quad D_G(x,y)\le t.
\]

Set \(t=19/10=1.9\). The edge \(ab\) is present because

\[
D_G(a,b)=\frac{3\sqrt{10}}5\approx1.897366596<1.9.
\]

Two other pairs are absent.

First, \(u\) has only the graph neighbors \(a,v\), so

\[
D_G(u,b)=
\min\left\{
\frac{\sqrt2+3\sqrt{10}}5,
\frac65+\frac{2}{\sqrt5}
\right\}
\approx2.094427191>1.9.
\tag{1}
\]

Both candidate paths have cost greater than \(2\): use \(\sqrt2>1\), \(\sqrt{10}>3\), and \(\sqrt5<5/2\). Thus the missing edge in (1) follows from exact inequalities.

Second, \(q\) has only the graph neighbors \(b,z\). Since \(ab\) and \(az\) are graph edges,

\[
D_G(a,q)=
\min\left\{\frac{4\sqrt{10}}5,\,2\sqrt2\right\}
=\frac{4\sqrt{10}}5
\approx2.529822128>1.9.
\tag{2}
\]

Now consider the two possible orientations of the existing Rips edge \(ab\):

- If \(a\to b\), the intermediate vertex \(u\) forces \(u\to b\) by (C), but \(ub\) is absent by (1).
- If \(b\to a\), the intermediate vertex \(q\) on \(b\prec q\prec z\prec a\) forces \(q\to a\), but \(qa\) is absent by (2).

Neither direction is possible. Therefore **\(H_{19/10}\) admits no cyclic orientation in the curve order**.

The obstruction holds throughout the interval

\[
\boxed{\frac{3\sqrt{10}}5\le t\le2},
\]

because \(ab\) is present and both obstructing pairs still have distance strictly greater than \(2\).

## 5. The obstruction occurs while the loop survives

At \(t=19/10\), direct shortest-path calculation gives

\[
E(H_t)=E(G)\cup\{uz,vq\}.
\]

Equivalently, the only new edges beyond \(G\) are \(uz\) and \(vq\), with distances

\[
D_G(u,z)=\frac{6\sqrt2}{5}\approx1.697056275,
\qquad
D_G(v,q)=\frac{2}{\sqrt5}+\frac{\sqrt{10}}5
\approx1.526882723.
\]

The four triangles are

\[
\{a,u,v\},\quad\{a,u,z\},\quad\{a,v,b\},\quad\{v,b,q\}.
\]

Their boundaries are independent over every field: respectively, the edges \(uv,uz,ab,vq\) each occur in only one triangle. The graph is connected and has ten edges and six vertices. Hence

\[
\dim H_1(\operatorname{VR}(D_G;19/10))
=10-6+1-4=1
\]

over any field. Thus the example has a surviving loop at the noncyclic threshold.

## 6. Normalization and scope

To use the paper's unit-circumference convention, divide every coordinate and distance by \(2\pi\). The threshold becomes

\[
t_{\mathrm{normalized}}=\frac{19}{20\pi}.
\]

The normalized cyclic arc coordinates are

\[
0,\quad\frac{\arctan(1/7)}\pi,\quad\frac14,
\quad\frac{\arctan3}\pi,\quad\frac12,\quad\frac34.
\]

The largest consecutive arc gap is \(g=1/4<1/3\), so the sample also lies in the paper's sufficient positive-reference regime. Uniform rescaling preserves arc support and the threshold graph.

The repeated chord lengths in these convenient coordinates are not essential. All support inequalities involving interior vertices and consecutive steps are strict, and the three threshold comparisons used in the obstruction are strict. Small generic angular perturbations therefore preserve the example and can make all unordered-pair chord lengths distinct.

The mechanism is explicit: arc support for \(ab\) supplies a fan from \(a\), but it need not make \(u\) close to \(b\). Cyclicity would require that additional pair in the forward direction. The missing pair \(aq\) blocks the reverse direction.

Consequently, arc support alone does not imply cyclicity of the graph-metric Rips threshold graphs in the curve order. This construction uses selected chords; the separate prefix results for radius, OR and mutual circle graphs retain their stated assumptions.

## 7. Reproducible verification

The following standalone Python code needs only NumPy. It checks all eight arc-support certificates, both missing pairs, the full threshold graph and triangle list, all \(2^{10}=1024\) orientations in the given curve order, and \(\beta_1=1\) over \(\mathbb F_2\). The proofs above establish the exact claims without relying on floating-point persistence output.

```python
import itertools
import numpy as np

# Indices 0,...,5 follow the geometric curve order a,u,v,b,q,z.
names = ['a', 'u', 'v', 'b', 'q', 'z']
points = np.array([
    [1., 0.], [24/25, 7/25], [0., 1.],
    [-4/5, 3/5], [-1., 0.], [0., -1.],
])
n = len(points)
ring = [(i, (i + 1) % n) for i in range(n)]
edges = ring + [(0, 2), (0, 3)]
chords = np.linalg.norm(points[:, None, :] - points[None, :, :], axis=2)
D = np.full((n, n), np.inf)
np.fill_diagonal(D, 0.)
for a, b in edges:
    D[a, b] = D[b, a] = chords[a, b]
for k in range(n):  # Floyd-Warshall, independent of the release's Dijkstra path.
    D = np.minimum(D, D[:, [k]] + D[[k], :])


def supporting_arc(a, b):
    cost = chords[a, b]
    for direction in (1, -1):
        steps = ((b - a) * direction) % n
        arc = (a + direction * np.arange(steps + 1)) % n
        if np.any(chords[arc[:-1], arc[1:]] > cost + 1e-12):
            continue
        for anchor in (a, b):
            if np.all(D[anchor, arc] <= cost + 1e-12):
                return arc.tolist(), anchor
    return None


assert all(supporting_arc(a, b) is not None for a, b in edges)
t = 19/10
H = [(a, b) for a, b in itertools.combinations(range(n), 2) if D[a, b] <= t]
H_set = set(H)
expected = {tuple(sorted(e)) for e in edges} | {(1, 5), (2, 4)}
assert H_set == expected
assert D[0, 3] < t < min(D[1, 3], D[0, 4])


def cyclic_in_given_order(directed):
    for a, b in directed:
        for step in range(1, (b - a) % n):
            middle = (a + step) % n
            if (a, middle) not in directed or (middle, b) not in directed:
                return False
    return True


valid_orientations = 0
for choices in itertools.product((0, 1), repeat=len(H)):
    directed = {(b, a) if reverse else (a, b)
                for (a, b), reverse in zip(H, choices)}
    valid_orientations += cyclic_in_given_order(directed)
assert valid_orientations == 0

triangles = [tri for tri in itertools.combinations(range(n), 3)
             if all(pair in H_set for pair in itertools.combinations(tri, 2))]
assert triangles == [(0, 1, 2), (0, 1, 5), (0, 2, 3), (2, 3, 4)]
edge_index = {e: i for i, e in enumerate(H)}
basis = {}
for triangle in triangles:
    column = sum(1 << edge_index[pair]
                 for pair in itertools.combinations(triangle, 2))
    while column:
        pivot = column.bit_length() - 1
        if pivot not in basis:
            basis[pivot] = column
            break
        column ^= basis[pivot]
beta1 = len(H) - n + 1 - len(basis)  # H is connected: it contains the ring.
assert beta1 == 1

print('Arc-supported edges:', len(edges), '/', len(edges))
print('Threshold:', t)
print('D(a,b), D(u,b), D(a,q):', D[0, 3], D[1, 3], D[0, 4])
print('Rips edges:', len(H))
print('Cyclic orientations in curve order:', valid_orientations)
print('H1 dimension over F2:', beta1)
```

Verified output:

```text
Arc-supported edges: 8 / 8
Threshold: 1.9
D(a,b), D(u,b), D(a,q): 1.8973665961010275 2.094427190999916 2.5298221281347035
Rips edges: 10
Cyclic orientations in curve order: 0
H1 dimension over F2: 1
```

## 8. No common cyclic order throughout the filtration

At \(t=6/5\), the threshold graph is exactly the path \(a-u-v-b-q\) and the isolated vertex \(z\). It is triangle-free. If an oriented edge spans an interior vertex in a cyclic order, condition (C) forces two additional edges, producing a triangle. Thus each path edge must join consecutive vertices in any admissible order. The path forces the curve order \((a,u,v,b,q,z)\) or its reverse, up to rotation.

At \(t=19/10\), neither orientation of \(ab\) is possible in that order, as Section 4 proves. Existence of a cyclic orientation is unchanged by reversing the order. Therefore **no single cyclic vertex order makes all threshold graphs cyclic**.

This is a filtration obstruction. The graph at \(t=19/10\) is cyclic in the different order \((a,u,z,q,b,v)\); it is false that this individual graph is noncyclic in every order. Indeed each distinct threshold graph has an admissible order separately.

The permanent checker `code/verify_separation.py` verifies all 120 orders (rotations identified; reversals retained) at all 14 distinct threshold graphs, including the empty graph. It solves the orientation implications and directly checks every returned orientation. At both proof graphs it also checks all directions by exhaustive enumeration, independently of that solver.

Metric thresholds are exact linear combinations in the basis \((1,\sqrt2,\sqrt5,\sqrt{10})/5\). Equal coefficient vectors identify tied distances; outward integer square-root bounds at 80 decimal places certify all comparisons. The certificate `results/separation_verification.json` lists all admissible orders and orientations at each threshold, together with their empty intersection.

The curve-order obstruction extends to the exact half-open interval

\[
\frac{3\sqrt{10}}5\le t<\frac65+\frac2{\sqrt5}.
\]

Its upper endpoint is \(D_G(u,b)\), so it is excluded. The two proof thresholds suffice for the common-order conclusion.
