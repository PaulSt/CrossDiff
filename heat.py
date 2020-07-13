import sys
from netgen.csg import unit_cube
from ngsolve.TensorProductTools import SegMesh
from netgen.geom2d import unit_square
from ngsolve.solvers import *
from prodmesh import TensorProdMesh,CartSquare
from ngsolve import *
from math import pi

def SolveHeat(order,eps,mesh,bndcond):
    V = H1(mesh, order=order)

    u0 = GridFunction(V)
    u0.Set(bndcond, BND, definedon=mesh.Boundaries('bottom'))

    D=mesh.dim
    w,v = V.TnT()
    u = exp(w)/(1+exp(w))
    dtv= CoefficientFunction(grad(v)[D-1])
    gv = CoefficientFunction(tuple([grad(v)[i] for i in  range(D-1)]))
    gw = CoefficientFunction(tuple([grad(w)[i] for i in  range(D-1)]))

    a = BilinearForm(V)
    a += SymbolicBFI(-u*dtv*tau + (u-u**2)*InnerProduct(gv,gw)+eps*w*v)
    a += SymbolicBFI(u*v*tau, BND, definedon=mesh.Boundaries('top'))
    a += SymbolicBFI(-u0*v*tau, BND, definedon=mesh.Boundaries('bottom'))

    gfu = GridFunction(V)
    with TaskManager():
        [success, steps]=Newton(a,gfu)
    return gfu

if __name__ == "__main__":
    eps=0
    D=2
    order = 4
    quads = True
    ratio = 2
    tau = 7
    shift=1.5
    scale=2.5*shift
    freq=1

    if D is 3:
        mesh = Mesh (unit_cube.GenerateMesh(maxh=0.5, quad_dominated=quads))
        bndcond = exp(-2*(freq*pi)**2*z/tau)*cos(freq*pi*x)*cos(freq*pi*y)
    else:
        mesh = Mesh (unit_square.GenerateMesh(maxh=0.5, quad_dominated=quads))
        bndcond = exp(-(freq*pi)**2*y/tau)*cos(freq*pi*x)
    bndcondscale = (bndcond+shift)/scale

    Maxh=[1/2**i for i in range(1,7)]
    errorold=0
    for i,maxh in enumerate(Maxh):
        if quads and D is 3:
            meshx = Mesh(unit_square.GenerateMesh(maxh=maxh))
            mesht = Mesh(SegMesh(round(1/(maxh)),0,1,periodic=False) )
            mesh = Mesh(TensorProdMesh(meshx,mesht))
        if quads and D is 2:
            # mesh = CartSquare(round(1/(maxh*2)),round(1/(maxh)))
            mesh = CartSquare(2**i,2**(i+1))

        gfu = SolveHeat(order,eps,mesh,bndcondscale)
        error = sqrt(Integrate((exp(gfu)/(1+exp(gfu))-bndcondscale)**2,mesh))
        rate=log(errorold/error)/log(2)
        errorold=error
        if not quads:
            mesh.Refine()
        print("===== order: ",order,' error: ', error, ' rate: ', rate)
        input()

