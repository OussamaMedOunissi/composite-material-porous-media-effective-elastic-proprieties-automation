import gmsh
import math
import sys
import numpy as nm
from shapely import geometry
import mesh_con
import matplotlib.pyplot as plt

#----- helpful functions ----------

def n_diff_area(s1,n,var):
    s = []
    while len(s) < n:
        b = min(20, n - len(s))
        t = nm.random.randint(1, b)
        while n - (t + len(s)) == 2:
            t = nm.random.randint(1, b)
        k = abs(nm.random.normal(0, var))
        while k >= t * 0.8:
            k = abs(nm.random.normal(0, var))
        s.append(s1 + k * s1)
        bh = min(0.8, k)
        bl = max(0, k - (0.8 * (t - 1)))
        h = nm.random.uniform(bl, bh)
        for j in range(t - 1):
            s.append(s1 - h * s1)
            k = k - h
            bh = min(0.8, k)
            bl = max(0, k - (0.8 * (t - j - 2)))
            h = nm.random.uniform(bh, bl)
        s.append(s1 - k * s1)
    return s

def n_diff_f(n,fmax,varf):
    f = []
    for i in range(n):
        k = abs(nm.random.normal(0, varf))
        while k > fmax - 1:
            k = abs(nm.random.normal(0, varf))
        f.append(k + 1)
    return f

def ellipse_inter(ellipses, n=100):
    t = nm.linspace(0, 2*nm.pi, n, endpoint=False)
    st = nm.sin(t)
    ct = nm.cos(t)
    result = [[],[]]
    for x0, y0, a, b, angle in ellipses:
        result[1].append([x0,y0])
        a *= 1
        b *= 1
        sa = nm.sin(angle)
        ca = nm.cos(angle)
        p = nm.empty((n, 2))
        p[:, 0] = x0 + a * ca * ct - b * sa * st
        p[:, 1] = y0 + a * sa * ct + b * ca * st
        p = nm.append(p, [p[0]], axis=0)
        result[0].append(p)
    e1, e2 = result[0]
    c1, c2 = result[1]
    ea = geometry.LineString(e1)
    eb = geometry.LineString(e2)
    p1 = geometry.Polygon(ea)
    p2 = geometry.Polygon(eb)
    o1 = geometry.Point(c1)
    o2 = geometry.Point(c2)
    t1 = ea.intersects(eb)
    t2 = p1.contains(o2)
    t3 = p2.contains(o1)
    return t1 or t2 or t3

def rect_inter(rects):
    result = [[],[]]
    for c, r, an in rects:
        result[1].append(c)
        l = (r[0] ** 2 + r[1] ** 2) ** (0.5)
        ang = []
        ang.append(an + nm.arctan(r[1] / r[0]))
        ang.append(an + (nm.pi - nm.arctan(r[1] / r[0])))
        ang.append(an - (nm.pi - nm.arctan(r[1] / r[0])))
        ang.append(an - nm.arctan(r[1] / r[0]))
        p = []
        for angl in ang:
            p.append([c[0] + l * nm.cos(angl),c[1] + l * nm.sin(angl)])
        p.append(p[0])
        pt = nm.array(p)
        result[0].append(pt)
    e1, e2 = result[0]
    c1, c2 = result[1]
    ea = geometry.LineString(e1)
    eb = geometry.LineString(e2)
    p1 = geometry.Polygon(ea)
    p2 = geometry.Polygon(eb)
    o1 = geometry.Point(c1)
    o2 = geometry.Point(c2)
    t1 = ea.intersects(eb)
    t2 = p1.contains(o2)
    t3 = p2.contains(o1)
    return t1 or t2 or t3

def trian_inter(trians):
    result = [[],[]]
    for c, r, an in trians:
        result[1].append(c)
        p = []
        ang = an + nm.pi/2
        p.append([c[0], c[1]])
        p.append([c[0] + r[0] * nm.cos(an), c[1] + r[0] * nm.sin(an)])
        p.append([c[0] + r[1] * nm.cos(ang), c[1] + r[1] * nm.sin(ang)])
        p.append(p[0])
        pt = nm.array(p)
        result[0].append(pt)
    e1, e2 = result[0]
    c1, c2 = result[1]
    ea = geometry.LineString(e1)
    eb = geometry.LineString(e2)
    p1 = geometry.Polygon(ea)
    p2 = geometry.Polygon(eb)
    o1 = geometry.Point(c1)
    o2 = geometry.Point(c2)
    t1 = ea.intersects(eb)
    t2 = p1.contains(o2)
    t3 = p2.contains(o1)
    return t1 or t2 or t3

def ellipse_ax(ar,f):
    prod = ar / nm.pi
    return [(prod * f) ** (0.5), (prod / f) ** (0.5)]

def rect_ax(ar,f):
    prod = ar
    return [((prod * f) ** (0.5))/2, ((prod / f) ** (0.5))/2]

def trian_ax(ar,f):
    prod = ar * 2
    return [((prod * f) ** (0.5)), ((prod / f) ** (0.5))]

#----- geometry of basic elements ----------

def create_square(p,lc):
    gmsh.model.geo.addPoint(p, p, 0, lc, 1)
    gmsh.model.geo.addPoint(-p, p, 0, lc, 2)
    gmsh.model.geo.addPoint(-p, -p, 0, lc, 3)
    gmsh.model.geo.addPoint(p, -p, 0, lc, 4)
    gmsh.model.geo.addLine(1, 2, 1)
    gmsh.model.geo.addLine(2, 3, 2)
    gmsh.model.geo.addLine(3, 4, 3)
    gmsh.model.geo.addLine(4, 1, 4)
    gmsh.model.geo.addCurveLoop([1, 2, 3, 4], 1)

def create_circle(cid,c,r,lc):
    gmsh.model.geo.addPoint(c[0], c[1], 0, lc, cid + 1)
    gmsh.model.geo.addPoint(c[0] + r, c[1], 0, lc, cid + 2)
    gmsh.model.geo.addPoint(c[0], c[1] + r, 0, lc, cid + 3)
    gmsh.model.geo.addPoint(c[0] - r, c[1], 0, lc, cid + 4)
    gmsh.model.geo.addPoint(c[0], c[1] - r, 0, lc, cid + 5)
    gmsh.model.geo.addCircleArc(cid + 2, cid + 1, cid + 3, cid + 1)
    gmsh.model.geo.addCircleArc(cid + 3, cid + 1, cid + 4, cid + 2)
    gmsh.model.geo.addCircleArc(cid + 4, cid + 1, cid + 5, cid + 3)
    gmsh.model.geo.addCircleArc(cid + 5, cid + 1, cid + 2, cid + 4)
    gmsh.model.geo.addCurveLoop([cid + 1, cid + 2, cid + 3, cid + 4], cid)
    gmsh.model.geo.addPlaneSurface([cid], cid)

def create_ellipse(cid,c,r,an,lc):
    gmsh.model.geo.addPoint(c[0], c[1], 0, lc, cid + 1)
    gmsh.model.geo.addPoint(c[0] + r[0]*nm.cos(an), c[1] + r[0]*nm.sin(an), 0, lc, cid + 2)
    gmsh.model.geo.addPoint(c[0] + r[1]*nm.cos((nm.pi/2) + an), c[1] + r[1]*nm.sin((nm.pi/2) + an), 0, lc, cid + 3)
    gmsh.model.geo.addPoint(c[0] - r[0]*nm.cos(an), c[1] -r[0]*nm.sin(an) , 0, lc, cid + 4)
    gmsh.model.geo.addPoint(c[0] - r[1]*nm.cos((nm.pi/2) + an), c[1] - r[1]*nm.sin((nm.pi/2) + an), 0, lc, cid + 5)
    gmsh.model.geo.addEllipseArc(cid + 2, cid + 1, cid + 2, cid + 3, cid + 1)
    gmsh.model.geo.addEllipseArc(cid + 3, cid + 1, cid + 2, cid + 4, cid + 2)
    gmsh.model.geo.addEllipseArc(cid + 4, cid + 1, cid + 2, cid + 5, cid + 3)
    gmsh.model.geo.addEllipseArc(cid + 5, cid + 1, cid + 2, cid + 2, cid + 4)
    gmsh.model.geo.addCurveLoop([cid + 1, (cid + 2), cid + 3, (cid + 4)], cid)
    gmsh.model.geo.addPlaneSurface([cid], cid)

def create_spline_circle(cid,c,r,lc):
    t = r * nm.sin(nm.pi / 4)
    gmsh.model.geo.addPoint(c[0] + r, c[1], 0, lc, cid + 2)
    gmsh.model.geo.addPoint(c[0] + t, c[1] + t, 0, lc, cid + 3)
    gmsh.model.geo.addPoint(c[0], c[1] + r, 0, lc, cid + 4)
    gmsh.model.geo.addPoint(c[0] - t, c[1] + t, 0, lc, cid + 5)
    gmsh.model.geo.addPoint(c[0] - r, c[1], 0, lc, cid + 6)
    gmsh.model.geo.addPoint(c[0] - t, c[1] - t, 0, lc, cid + 7)
    gmsh.model.geo.addPoint(c[0], c[1] - r, 0, lc, cid + 8)
    gmsh.model.geo.addPoint(c[0] + t, c[1] - t, 0, lc, cid + 9)
    gmsh.model.geo.addSpline([cid + 2, cid + 3, cid + 4, cid + 5, cid + 6, cid + 7, cid + 8, cid + 9, cid + 2], cid + 1)
    gmsh.model.geo.addCurveLoop([cid + 1], cid)
    gmsh.model.geo.addPlaneSurface([cid], cid)

def create_rectangle(cid,c,r,an,lc):
    l = (r[0]**2 + r[1]**2)**(0.5)
    ang1 = an + nm.arctan(r[1]/r[0])
    ang2 = an + (nm.pi - nm.arctan(r[1] / r[0]))
    ang3 = an - (nm.pi - nm.arctan(r[1] / r[0]))
    ang4 = an - nm.arctan(r[1] / r[0])
    gmsh.model.geo.addPoint(c[0] + l*nm.cos(ang1), c[1] + l*nm.sin(ang1), 0, lc, cid + 2)
    gmsh.model.geo.addPoint(c[0] + l*nm.cos(ang2), c[1] + l*nm.sin(ang2), 0, lc, cid + 3)
    gmsh.model.geo.addPoint(c[0] + l*nm.cos(ang3), c[1] + l*nm.sin(ang3), 0, lc, cid + 4)
    gmsh.model.geo.addPoint(c[0] + l*nm.cos(ang4), c[1] + l*nm.sin(ang4), 0, lc, cid + 5)
    gmsh.model.geo.addLine(cid + 2, cid + 3, cid + 1)
    gmsh.model.geo.addLine(cid + 3, cid + 4, cid + 2)
    gmsh.model.geo.addLine(cid + 4, cid + 5, cid + 3)
    gmsh.model.geo.addLine(cid + 5, cid + 2, cid + 4)
    gmsh.model.geo.addCurveLoop([cid + 1, (cid + 2), cid + 3, (cid + 4)], cid)
    gmsh.model.geo.addPlaneSurface([cid], cid)

def create_actriangle(cid,c,r,an,lc):
    ang = an + nm.pi/2
    gmsh.model.geo.addPoint(c[0], c[1], 0, lc, cid + 1)
    gmsh.model.geo.addPoint(c[0] + r[0]*nm.cos(an), c[1] + r[0]*nm.sin(an), 0, lc, cid + 2)
    gmsh.model.geo.addPoint(c[0] + r[1]*nm.cos(ang), c[1] + r[1]*nm.sin(ang), 0, lc, cid + 3)
    gmsh.model.geo.addLine(cid + 1, cid + 2, cid + 1)
    gmsh.model.geo.addLine(cid + 2, cid + 3, cid + 2)
    gmsh.model.geo.addLine(cid + 3, cid + 1, cid + 3)
    gmsh.model.geo.addCurveLoop([cid + 1, (cid + 2), cid + 3], cid)
    gmsh.model.geo.addPlaneSurface([cid], cid)

#----- parameters of different composite materials structure ----------

# ---- cercle ----

def n_cercle(n,p,lc):
    r = [((p * 4) / (nm.pi * n)) ** (0.5)] * n
    c = [[nm.random.uniform(-1 + r[0] + lc, 1 - r[0] - lc), nm.random.uniform(-1 + r[0] + lc, 1 - r[0] - lc)]]
    for i in range(n - 1):
        vi = [nm.random.uniform(-1 + r[i + 1] + lc, 1 - r[i + 1] - lc),
              nm.random.uniform(-1 + r[i + 1] + lc, 1 - r[i + 1] - lc)]
        inter = True
        while inter:
            f = True
            for v in range(len(c)):
                if ((c[v][0] - vi[0]) ** 2 + (c[v][1] - vi[1]) ** 2) ** (0.5) <= r[i + 1] + r[v]:
                    vi = [nm.random.uniform(-1 + r[i + 1] + lc, 1 - r[i + 1] - lc),
                          nm.random.uniform(-1 + r[i + 1] + lc, 1 - r[i + 1] - lc)]
                    f = False
                    break
            if f:
                inter = False
        c.append(vi)
    return [n,c,r]

def n_diff_cercle(n,p,lc,var):
    s1 = ((p * 4) / (n))
    s = n_diff_area(s1, n, var)
    r = []
    for e in s:
        r.append((e/nm.pi)**(0.5))
    c = [[nm.random.uniform(-1+r[0]+lc, 1-r[0]-lc), nm.random.uniform(-1+r[0]+lc, 1-r[0]-lc)]]
    for i in range(n-1):
        vi = [nm.random.uniform(-1+r[i+1]+lc, 1-r[i+1]-lc), nm.random.uniform(-1+r[i+1]+lc, 1-r[i+1]-lc)]
        inter = True
        while inter:
            f = True
            for v in range(len(c)):
                if ((c[v][0] - vi[0])**2 + (c[v][1] - vi[1])**2)**(0.5) <= r[i+1]+r[v]:
                    vi = [nm.random.uniform(-1+r[i+1]+lc, 1-r[i+1]-lc), nm.random.uniform(-1+r[i+1]+lc, 1-r[i+1]-lc)]
                    f = False
                    break
            if f:
                inter = False
        c.append(vi)
    return [n,c,r]

# ---- ellipse ----

def n_ellipse(n,p,ff,an,lc):
    ar = ((p * 4) / (n))
    r = ellipse_ax(ar,ff)
    ry = [r] * n
    ang = [an] * n
    c = [[nm.random.uniform(-1+r[0] + lc, 1-r[0]-lc), nm.random.uniform(-1+r[0]+lc, 1-r[0]-lc)]]
    for i in range(n-1):
        vi = [nm.random.uniform(-1+r[0]+lc, 1-r[0]-lc), nm.random.uniform(-1+r[0]+lc, 1-r[0]-lc)]
        inter = True
        while inter:
            f = True
            for v in c:
                if ellipse_inter([(*v,*r,an),(*vi,*r,an)]):
                    vi = [nm.random.uniform(-1+r[0]+lc, 1-r[0]-lc), nm.random.uniform(-1+r[0]+lc, 1-r[0]-lc)]
                    f = False
                    break
            if f:
                inter = False
        c.append(vi)
    return [n, c, ry, ang]

def n_difan_ellipse(n,p,ff,lc):
    ar = ((p * 4) / (n))
    r = ellipse_ax(ar,ff)
    ry = [r] * n
    an = [nm.random.uniform(0, nm.pi)]
    c = [[nm.random.uniform(-1+r[0]+lc, 1-r[0]-lc), nm.random.uniform(-1+r[0]+lc, 1-r[0]-lc)]]
    for i in range(n-1):
        ani = nm.random.uniform(0, nm.pi)
        vi = [nm.random.uniform(-1+r[0]+lc, 1-r[0]-lc), nm.random.uniform(-1+r[0]+lc, 1-r[0]-lc)]
        inter = True
        while inter:
            f = True
            k = 0
            for v in c:
                if ellipse_inter([(*v,*r,an[k]),(*vi,*r,ani)]):
                    ani = nm.random.uniform(0, nm.pi)
                    vi = [nm.random.uniform(-1+r[0]+lc, 1-r[0]-lc), nm.random.uniform(-1+r[0]+lc, 1-r[0]-lc)]
                    f = False
                    break
                k += 1
            if f:
                inter = False
        c.append(vi)
        an.append(ani)
    return [n, c, ry, an]

def n_difanf_ellipse(n,p,fmax,lc,varf):
    ar = ((p * 4) / (n))
    fn = n_diff_f(n, fmax, varf)
    r = []
    for i in range(n):
        r.append(ellipse_ax(ar, fn[i]))
    an = [nm.random.uniform(0, nm.pi)]
    c = [[nm.random.uniform(-1+r[0][0]+lc, 1-r[0][0]-lc), nm.random.uniform(-1+r[0][0]+lc, 1-r[0][0]-lc)]]
    for i in range(n-1):
        ani = nm.random.uniform(0, nm.pi)
        vi = [nm.random.uniform(-1+r[i+1][0]+lc, 1-r[i+1][0]-lc), nm.random.uniform(-1+r[i+1][0]+lc, 1-r[i+1][0]-lc)]
        inter = True
        while inter:
            f = True
            k = 0
            for v in c:
                if ellipse_inter([(*v,*r[k],an[k]),(*vi,*r[i+1],ani)]):
                    ani = nm.random.uniform(0, nm.pi)
                    vi = [nm.random.uniform(-1+r[i+1][0]+lc, 1-r[i+1][0]-lc), nm.random.uniform(-1+r[i+1][0]+lc, 1-r[i+1][0]-lc)]
                    f = False
                    break
                k += 1
            if f:
                inter = False
        c.append(vi)
        an.append(ani)
    return [n,c,r,an]

def n_difanr_ellipse(n,p,fmax,lc,varr):
    ar = ((p * 4) / (n))
    s = n_diff_area(ar, n, varr)
    r = []
    for i in range(n):
        r.append(ellipse_ax(s[i], fmax))
    an = [nm.random.uniform(0, nm.pi)]
    c = [[nm.random.uniform(-1+r[0][0]+lc, 1-r[0][0]-lc), nm.random.uniform(-1+r[0][0]+lc, 1-r[0][0]-lc)]]
    for i in range(n-1):
        ani = nm.random.uniform(0, nm.pi)
        vi = [nm.random.uniform(-1+r[i+1][0]+lc, 1-r[i+1][0]-lc), nm.random.uniform(-1+r[i+1][0]+lc, 1-r[i+1][0]-lc)]
        inter = True
        while inter:
            f = True
            k = 0
            for v in c:
                if ellipse_inter([(*v,*r[k],an[k]),(*vi,*r[i+1],ani)]):
                    ani = nm.random.uniform(0, nm.pi)
                    vi = [nm.random.uniform(-1+r[i+1][0]+lc, 1-r[i+1][0]-lc), nm.random.uniform(-1+r[i+1][0]+lc, 1-r[i+1][0]-lc)]
                    f = False
                    break
                k += 1
            if f:
                inter = False
        c.append(vi)
        an.append(ani)
    return [n,c,r,an]

def n_difanfr_ellipse(n,p,fmax,lc,varr,varf):
    ar = ((p * 4) / (n))
    fn = n_diff_f(n, fmax,varf)
    s = n_diff_area(ar, n, varr)
    r = []
    for i in range(n):
        r.append(ellipse_ax(s[i], fn[i]))
    an = [nm.random.uniform(0, nm.pi)]
    c = [[nm.random.uniform(-1+r[0][0]+lc,  1-r[0][0]-lc), nm.random.uniform(-1+r[0][0]+lc, 1-r[0][0]-lc)]]
    for i in range(n-1):
        ani = nm.random.uniform(0, nm.pi)
        vi = [nm.random.uniform(-1+r[i+1][0]+lc, 1-r[i+1][0]-lc), nm.random.uniform(-1+r[i+1][0]+lc, 1-r[i+1][0]-lc)]
        inter = True
        while inter:
            f = True
            k = 0
            for v in c:
                if ellipse_inter([(*v,*r[k],an[k]),(*vi,*r[i+1],ani)]):
                    ani = nm.random.uniform(0, nm.pi)
                    vi = [nm.random.uniform(-1+r[i+1][0]+lc, 1-r[i+1][0]-lc), nm.random.uniform(-1+r[i+1][0]+lc, 1-r[i+1][0]-lc)]
                    f = False
                    break
                k += 1
            if f:
                inter = False
        c.append(vi)
        an.append(ani)
    return [n,c,r,an]

# ---- rectangle ----

def n_rectangle(n,p,ff,an,lc):
    ar = ((p * 4) / (n))
    r = rect_ax(ar,ff)
    ry = [r] * n
    ang = [an] * n
    c = [[nm.random.uniform(-1+r[0] + lc, 1-r[0]-lc), nm.random.uniform(-1+r[0]+lc, 1-r[0]-lc)]]
    for i in range(n-1):
        vi = [nm.random.uniform(-1+r[0]+lc, 1-r[0]-lc), nm.random.uniform(-1+r[0]+lc, 1-r[0]-lc)]
        inter = True
        while inter:
            f = True
            for v in c:
                if rect_inter([(v,r,an),(vi,r,an)]):
                    vi = [nm.random.uniform(-1+r[0]+lc, 1-r[0]-lc), nm.random.uniform(-1+r[0]+lc, 1-r[0]-lc)]
                    f = False
                    break
            if f:
                inter = False
        c.append(vi)
    return [n, c, ry, ang]

def n_difan_rectangle(n,p,ff,lc):
    ar = ((p * 4) / (n))
    r = rect_ax(ar,ff)
    ry = [r] * n
    an = [nm.random.uniform(0, nm.pi)]
    c = [[nm.random.uniform(-1+r[0]+lc, 1-r[0]-lc), nm.random.uniform(-1+r[0]+lc, 1-r[0]-lc)]]
    for i in range(n-1):
        ani = nm.random.uniform(0, nm.pi)
        vi = [nm.random.uniform(-1+r[0]+lc, 1-r[0]-lc), nm.random.uniform(-1+r[0]+lc, 1-r[0]-lc)]
        inter = True
        while inter:
            f = True
            k = 0
            for v in c:
                if rect_inter([(v,r,an[k]),(vi,r,ani)]):
                    ani = nm.random.uniform(0, nm.pi)
                    vi = [nm.random.uniform(-1+r[0]+lc, 1-r[0]-lc), nm.random.uniform(-1+r[0]+lc, 1-r[0]-lc)]
                    f = False
                    break
                k += 1
            if f:
                inter = False
        c.append(vi)
        an.append(ani)
    return [n, c, ry, an]

def n_difanf_rectangle(n,p,fmax,lc,varf):
    ar = ((p * 4) / (n))
    fn = n_diff_f(n, fmax, varf)
    r = []
    for i in range(n):
        r.append(rect_ax(ar, fn[i]))
    an = [nm.random.uniform(0, nm.pi)]
    c = [[nm.random.uniform(-1+r[0][0]+lc, 1-r[0][0]-lc), nm.random.uniform(-1+r[0][0]+lc, 1-r[0][0]-lc)]]
    for i in range(n-1):
        ani = nm.random.uniform(0, nm.pi)
        vi = [nm.random.uniform(-1+r[i+1][0]+lc, 1-r[i+1][0]-lc), nm.random.uniform(-1+r[i+1][0]+lc, 1-r[i+1][0]-lc)]
        inter = True
        while inter:
            f = True
            k = 0
            for v in c:
                if rect_inter([(v,r[k],an[k]),(vi,r[i+1],ani)]):
                    ani = nm.random.uniform(0, nm.pi)
                    vi = [nm.random.uniform(-1+r[i+1][0]+lc, 1-r[i+1][0]-lc), nm.random.uniform(-1+r[i+1][0]+lc, 1-r[i+1][0]-lc)]
                    f = False
                    break
                k += 1
            if f:
                inter = False
        c.append(vi)
        an.append(ani)
    return [n,c,r,an]

def n_difanr_rectangle(n,p,fmax,lc,varr):
    ar = ((p * 4) / (n))
    s = n_diff_area(ar, n, varr)
    r = []
    for i in range(n):
        r.append(rect_ax(s[i], fmax))
    an = [nm.random.uniform(0, nm.pi)]
    c = [[nm.random.uniform(-1+r[0][0]+lc, 1-r[0][0]-lc), nm.random.uniform(-1+r[0][0]+lc, 1-r[0][0]-lc)]]
    for i in range(n-1):
        ani = nm.random.uniform(0, nm.pi)
        vi = [nm.random.uniform(-1+r[i+1][0]+lc, 1-r[i+1][0]-lc), nm.random.uniform(-1+r[i+1][0]+lc, 1-r[i+1][0]-lc)]
        inter = True
        while inter:
            f = True
            k = 0
            for v in c:
                if rect_inter([(v,r[k],an[k]),(vi,r[i+1],ani)]):
                    ani = nm.random.uniform(0, nm.pi)
                    vi = [nm.random.uniform(-1+r[i+1][0]+lc, 1-r[i+1][0]-lc), nm.random.uniform(-1+r[i+1][0]+lc, 1-r[i+1][0]-lc)]
                    f = False
                    break
                k += 1
            if f:
                inter = False
        c.append(vi)
        an.append(ani)
    return [n,c,r,an]

def n_difanfr_rectangle(n,p,fmax,lc,varr,varf):
    ar = ((p * 4) / (n))
    fn = n_diff_f(n, fmax,varf)
    s = n_diff_area(ar, n, varr)
    r = []
    for i in range(n):
        r.append(rect_ax(s[i], fn[i]))
    an = [nm.random.uniform(0, nm.pi)]
    c = [[nm.random.uniform(-1+r[0][0]+lc,  1-r[0][0]-lc), nm.random.uniform(-1+r[0][0]+lc, 1-r[0][0]-lc)]]
    for i in range(n-1):
        ani = nm.random.uniform(0, nm.pi)
        vi = [nm.random.uniform(-1+r[i+1][0]+lc, 1-r[i+1][0]-lc), nm.random.uniform(-1+r[i+1][0]+lc, 1-r[i+1][0]-lc)]
        inter = True
        while inter:
            f = True
            k = 0
            for v in c:
                if rect_inter([(v,r[k],an[k]),(vi,r[i+1],ani)]):
                    ani = nm.random.uniform(0, nm.pi)
                    vi = [nm.random.uniform(-1+r[i+1][0]+lc, 1-r[i+1][0]-lc), nm.random.uniform(-1+r[i+1][0]+lc, 1-r[i+1][0]-lc)]
                    f = False
                    break
                k += 1
            if f:
                inter = False
        c.append(vi)
        an.append(ani)
    return [n,c,r,an]

# ---- triangle ----

def n_triangle(n,p,ff,an,lc):
    ar = ((p * 4) / (n))
    r = trian_ax(ar,ff)
    ry = [r] * n
    ang = [an] * n
    c = [[nm.random.uniform(-1+r[0] + lc, 1-r[0]-lc), nm.random.uniform(-1+r[0]+lc, 1-r[0]-lc)]]
    for i in range(n-1):
        vi = [nm.random.uniform(-1+r[0]+lc, 1-r[0]-lc), nm.random.uniform(-1+r[0]+lc, 1-r[0]-lc)]
        inter = True
        while inter:
            f = True
            for v in c:
                if trian_inter([(v,r,an),(vi,r,an)]):
                    vi = [nm.random.uniform(-1+r[0]+lc, 1-r[0]-lc), nm.random.uniform(-1+r[0]+lc, 1-r[0]-lc)]
                    f = False
                    break
            if f:
                inter = False
        c.append(vi)
    return [n, c, ry, ang]

def n_difan_triangle(n,p,ff,lc):
    ar = ((p * 4) / (n))
    r = trian_ax(ar,ff)
    ry = [r] * n
    an = [nm.random.uniform(0, nm.pi*2)]
    c = [[nm.random.uniform(-1+r[0]+lc, 1-r[0]-lc), nm.random.uniform(-1+r[0]+lc, 1-r[0]-lc)]]
    for i in range(n-1):
        ani = nm.random.uniform(0, nm.pi*2)
        vi = [nm.random.uniform(-1+r[0]+lc, 1-r[0]-lc), nm.random.uniform(-1+r[0]+lc, 1-r[0]-lc)]
        inter = True
        while inter:
            f = True
            k = 0
            for v in c:
                if trian_inter([(v,r,an[k]),(vi,r,ani)]):
                    ani = nm.random.uniform(0, nm.pi*2)
                    vi = [nm.random.uniform(-1+r[0]+lc, 1-r[0]-lc), nm.random.uniform(-1+r[0]+lc, 1-r[0]-lc)]
                    f = False
                    break
                k += 1
            if f:
                inter = False
        c.append(vi)
        an.append(ani)
    return [n, c, ry, an]

#----- geometry of n structure ----------

def mesh_n_cercle(n,c,r,lc):
    civ = []
    for j in range(n):
        create_circle((j + 1) * 100, c[j], r[j], lc)
        civ.append((j + 1) * 100)
    return civ

def mesh_n_ellipse(n,c,r,an,lc):
    civ = []
    for j in range(n):
        create_ellipse((j+1)*100, c[j], r[j], an[j], lc)
        civ.append((j+1)*100)
    return civ

def mesh_n_rectangle(n,c,r,an,lc):
    civ = []
    for j in range(n):
        create_rectangle((j+1)*100, c[j], r[j], an[j], lc)
        civ.append((j+1)*100)
    return civ

def mesh_n_triangle(n,c,r,an,lc):
    civ = []
    for j in range(n):
        create_actriangle((j+1)*100, c[j], r[j], an[j], lc)
        civ.append((j+1)*100)
    return civ

#----- mesh the different structure ----------

# ---- cercle ----

def create_mesh_n_cercle(names,lc,e,p,n):
    data = n_cercle(n, p, max(lc))
    i = 0
    for name in names:
        gmsh.initialize()
        gmsh.model.add(name)
        create_square(e, lc[i])
        civ = mesh_n_cercle(*data, lc[i])
        gmsh.model.geo.addPlaneSurface([1, *civ], 1)
        gmsh.model.geo.synchronize()
        gmsh.model.addPhysicalGroup(2, [1], 1)
        gmsh.model.addPhysicalGroup(2, civ, 2)
        gmsh.model.mesh.generate(2)
        gmsh.write(name + ".mesh")
        mesh_con.convert_2d(name + ".mesh", name + ".h5")
        gmsh.fltk.run()
        gmsh.finalize()
        i += 1

def create_mesh_diff_cercle(names,lc,e,p,n,var):
    data = n_diff_cercle(n, p, max(lc),var)
    i = 0
    for name in names:
        gmsh.initialize()
        gmsh.model.add(name)
        create_square(e, lc[i])
        civ = mesh_n_cercle(*data, lc[i])
        gmsh.model.geo.addPlaneSurface([1, *civ], 1)
        gmsh.model.geo.synchronize()
        gmsh.model.addPhysicalGroup(2, [1], 1)
        gmsh.model.addPhysicalGroup(2, civ, 2)
        gmsh.model.mesh.generate(2)
        gmsh.write(name + ".mesh")
        mesh_con.convert_2d(name + ".mesh", name + ".h5")
        gmsh.fltk.run()
        gmsh.finalize()
        i += 1

# ---- ellipse ----

def create_mesh_n_ellipse(names,lc,e,n,p,ff,an):
    data = n_ellipse(n,p,ff,an,max(lc))
    i = 0
    for name in names:
        gmsh.initialize()
        gmsh.model.add(name)
        create_square(e, lc[i])
        civ = mesh_n_ellipse(*data, lc[i])
        gmsh.model.geo.addPlaneSurface([1, *civ], 1)
        gmsh.model.geo.synchronize()
        gmsh.model.addPhysicalGroup(2, [1], 1)
        gmsh.model.addPhysicalGroup(2, civ, 2)
        gmsh.model.mesh.generate(2)
        gmsh.write(name + ".mesh")
        mesh_con.convert_2d(name + ".mesh", name + ".h5")
        gmsh.fltk.run()
        gmsh.finalize()
        i += 1

def create_mesh_n_difan_ellipse(names,lc,e,n,p,ff):
    data = n_difan_ellipse(n,p,ff,max(lc))
    i = 0
    for name in names:
        gmsh.initialize()
        gmsh.model.add(name)
        create_square(e, lc[i])
        civ = mesh_n_ellipse(*data, lc[i])
        gmsh.model.geo.addPlaneSurface([1, *civ], 1)
        gmsh.model.geo.synchronize()
        gmsh.model.addPhysicalGroup(2, [1], 1)
        gmsh.model.addPhysicalGroup(2, civ, 2)
        gmsh.model.mesh.generate(2)
        gmsh.write(name + ".mesh")
        mesh_con.convert_2d(name + ".mesh", name + ".h5")
        gmsh.fltk.run()
        gmsh.finalize()
        i += 1

def create_mesh_n_difanf_ellipse(names,lc,e,n,p,fmax,varf):
    data = n_difanf_ellipse(n,p,fmax,max(lc),varf)
    i = 0
    for name in names:
        gmsh.initialize()
        gmsh.model.add(name)
        create_square(e, lc[i])
        civ = mesh_n_ellipse(*data, lc[i])
        gmsh.model.geo.addPlaneSurface([1, *civ], 1)
        gmsh.model.geo.synchronize()
        gmsh.model.addPhysicalGroup(2, [1], 1)
        gmsh.model.addPhysicalGroup(2, civ, 2)
        gmsh.model.mesh.generate(2)
        gmsh.write(name + ".mesh")
        mesh_con.convert_2d(name + ".mesh", name + ".h5")
        gmsh.fltk.run()
        gmsh.finalize()
        i += 1

def create_mesh_n_difanr_ellipse(names,lc,e,n,p,fmax,varr):
    data = n_difanr_ellipse(n,p,fmax,max(lc),varr)
    i = 0
    for name in names:
        gmsh.initialize()
        gmsh.model.add(name)
        create_square(e, lc[i])
        civ = mesh_n_ellipse(*data, lc[i])
        gmsh.model.geo.addPlaneSurface([1, *civ], 1)
        gmsh.model.geo.synchronize()
        gmsh.model.addPhysicalGroup(2, [1], 1)
        gmsh.model.addPhysicalGroup(2, civ, 2)
        gmsh.model.mesh.generate(2)
        gmsh.write(name + ".mesh")
        mesh_con.convert_2d(name + ".mesh", name + ".h5")
        gmsh.fltk.run()
        gmsh.finalize()
        i += 1

def create_mesh_n_difanfr_ellipse(names,lc,e,n,p,fmax,varr,varf):
    data = n_difanfr_ellipse(n,p,fmax,max(lc),varr,varf)
    i = 0
    for name in names:
        gmsh.initialize()
        gmsh.model.add(name)
        create_square(e, lc[i])
        civ = mesh_n_ellipse(*data, lc[i])
        gmsh.model.geo.addPlaneSurface([1, *civ], 1)
        gmsh.model.geo.synchronize()
        gmsh.model.addPhysicalGroup(2, [1], 1)
        gmsh.model.addPhysicalGroup(2, civ, 2)
        gmsh.model.mesh.generate(2)
        gmsh.write(name + ".mesh")
        mesh_con.convert_2d(name + ".mesh", name + ".h5")
        gmsh.fltk.run()
        gmsh.finalize()
        i += 1

# ---- rectangle ----

def create_mesh_n_rectangle(names,lc,e,n,p,ff,an):
    data = n_rectangle(n,p,ff,an,max(lc))
    i = 0
    for name in names:
        gmsh.initialize()
        gmsh.model.add(name)
        create_square(e, lc[i])
        civ = mesh_n_rectangle(*data, lc[i])
        gmsh.model.geo.addPlaneSurface([1, *civ], 1)
        gmsh.model.geo.synchronize()
        gmsh.model.addPhysicalGroup(2, [1], 1)
        gmsh.model.addPhysicalGroup(2, civ, 2)
        gmsh.model.mesh.generate(2)
        gmsh.write(name + ".mesh")
        mesh_con.convert_2d(name + ".mesh", name + ".h5")
        gmsh.fltk.run()
        gmsh.finalize()
        i += 1

def create_mesh_n_difan_rectangle(names,lc,e,n,p,ff):
    data = n_difan_rectangle(n,p,ff,max(lc))
    i = 0
    for name in names:
        gmsh.initialize()
        gmsh.model.add(name)
        create_square(e, lc[i])
        civ = mesh_n_rectangle(*data, lc[i])
        gmsh.model.geo.addPlaneSurface([1, *civ], 1)
        gmsh.model.geo.synchronize()
        gmsh.model.addPhysicalGroup(2, [1], 1)
        gmsh.model.addPhysicalGroup(2, civ, 2)
        gmsh.model.mesh.generate(2)
        gmsh.write(name + ".mesh")
        mesh_con.convert_2d(name + ".mesh", name + ".h5")
        gmsh.fltk.run()
        gmsh.finalize()
        i += 1

def create_mesh_n_difanf_rectangle(names,lc,e,n,p,fmax,varf):
    data = n_difanf_rectangle(n,p,fmax,max(lc),varf)
    i = 0
    for name in names:
        gmsh.initialize()
        gmsh.model.add(name)
        create_square(e, lc[i])
        civ = mesh_n_rectangle(*data, lc[i])
        gmsh.model.geo.addPlaneSurface([1, *civ], 1)
        gmsh.model.geo.synchronize()
        gmsh.model.addPhysicalGroup(2, [1], 1)
        gmsh.model.addPhysicalGroup(2, civ, 2)
        gmsh.model.mesh.generate(2)
        gmsh.write(name + ".mesh")
        mesh_con.convert_2d(name + ".mesh", name + ".h5")
        gmsh.fltk.run()
        gmsh.finalize()
        i += 1

def create_mesh_n_difanr_rectangle(names,lc,e,n,p,fmax,varr):
    data = n_difanr_rectangle(n,p,fmax,max(lc),varr)
    i = 0
    for name in names:
        gmsh.initialize()
        gmsh.model.add(name)
        create_square(e, lc[i])
        civ = mesh_n_rectangle(*data, lc[i])
        gmsh.model.geo.addPlaneSurface([1, *civ], 1)
        gmsh.model.geo.synchronize()
        gmsh.model.addPhysicalGroup(2, [1], 1)
        gmsh.model.addPhysicalGroup(2, civ, 2)
        gmsh.model.mesh.generate(2)
        gmsh.write(name + ".mesh")
        mesh_con.convert_2d(name + ".mesh", name + ".h5")
        gmsh.fltk.run()
        gmsh.finalize()
        i += 1

def create_mesh_n_difanfr_rectangle(names,lc,e,n,p,fmax,varr,varf):
    data = n_difanfr_rectangle(n,p,fmax,max(lc),varr,varf)
    i = 0
    for name in names:
        gmsh.initialize()
        gmsh.model.add(name)
        create_square(e, lc[i])
        civ = mesh_n_rectangle(*data, lc[i])
        gmsh.model.geo.addPlaneSurface([1, *civ], 1)
        gmsh.model.geo.synchronize()
        gmsh.model.addPhysicalGroup(2, [1], 1)
        gmsh.model.addPhysicalGroup(2, civ, 2)
        gmsh.model.mesh.generate(2)
        gmsh.write(name + ".mesh")
        mesh_con.convert_2d(name + ".mesh", name + ".h5")
        gmsh.fltk.run()
        gmsh.finalize()
        i += 1

# ---- triangle ----

def create_mesh_n_triangle(names,lc,e,n,p,ff):
    data = n_triangle(n,p,ff,an,max(lc))
    i = 0
    for name in names:
        gmsh.initialize()
        gmsh.model.add(name)
        create_square(e, lc[i])
        civ = mesh_n_triangle(*data, lc[i])
        gmsh.model.geo.addPlaneSurface([1, *civ], 1)
        gmsh.model.geo.synchronize()
        gmsh.model.addPhysicalGroup(2, [1], 1)
        gmsh.model.addPhysicalGroup(2, civ, 2)
        gmsh.model.mesh.generate(2)
        gmsh.write(name + ".mesh")
        mesh_con.convert_2d(name + ".mesh", name + ".h5")
        gmsh.fltk.run()
        gmsh.finalize()
        i += 1

def create_mesh_n_difan_triangle(names,lc,e,n,p,ff):
    data = n_difan_triangle(n,p,ff,max(lc))
    i = 0
    for name in names:
        gmsh.initialize()
        gmsh.model.add(name)
        create_square(e, lc[i])
        civ = mesh_n_triangle(*data, lc[i])
        gmsh.model.geo.addPlaneSurface([1, *civ], 1)
        gmsh.model.geo.synchronize()
        gmsh.model.addPhysicalGroup(2, [1], 1)
        gmsh.model.addPhysicalGroup(2, civ, 2)
        gmsh.model.mesh.generate(2)
        gmsh.write(name + ".mesh")
        mesh_con.convert_2d(name + ".mesh", name + ".h5")
        gmsh.fltk.run()
        gmsh.finalize()
        i += 1

def create_name(type,lc,ni):
    names =[]
    for d in lc:
        names.append(type + str(d) + 'n' + str(ni))
    return names

# application
l = 2
e = 1
lc = [2e-2,1.5e-2]
p = [0.2]*l
n = [50,100,200,300,400,500,600]
fmax = [3]*l
ff = [1]*l
var = [1] * l
varr = [1] * l
varf = [1] * l
types = ["Input/cercle", "Input/ellipse", "Input/rectangle"]

for i in range(l):
    names = create_name(types[0],lc,i)
    create_mesh_diff_cercle(names,lc,e,p[i],n[i],var[i])
    names = create_name(types[1],lc,i)
    create_mesh_n_difanfr_ellipse(names,lc,e,n[i],p[i],fmax[i],varr[i],varf[i])
    names = create_name(types[2],lc,i)
    create_mesh_n_difan_rectangle(names,lc,e,n[i],p[i],ff[i])

name = "cs"
an = nm.pi/4
cid = 100
c = [0,0]
r = [0.5,0.5]

#create_mesh_n_cercle(names,lc,e,p,n)
#create_mesh_diff_cercle(names,lc,e,p,n,var)

#create_mesh_n_ellipse(names,lc,e,n,p,ff,an)
#create_mesh_n_difan_ellipse(names,lc,e,n,p,ff)
#create_mesh_n_difanf_ellipse(names,lc,e,n,p,fmax,varf)
#create_mesh_n_difanr_ellipse(names,lc,e,n,p,fmax,varr)
#create_mesh_n_difanfr_ellipse(names,lc,e,n,p,fmax,varr,varf)

#create_mesh_n_rectangle(names,lc,e,n,p,ff,an)
#create_mesh_n_difan_rectangle(names,lc,e,n,p,ff)
#create_mesh_n_difanf_rectangle(names,lc,e,n,p,fmax,varf)
#create_mesh_n_difanr_rectangle(names,lc,e,n,p,fmax,varr)
#create_mesh_n_difanfr_rectangle(names,lc,e,n,p,fmax,varf,varf)

#create_mesh_n_triangle(names,lc,e,n,p,ff)
#create_mesh_n_difan_triangle(names,lc,e,n,p,ff)