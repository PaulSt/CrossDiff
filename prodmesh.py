import netgen.meshing as ngm
from ngsolve import *
from netgen.geom2d import unit_square
from ngsolve.TensorProductTools import *
from ngsolve import Mesh
import netgen.geom2d as ngeom2d


def ProdMesh(ngmeshbase,t_steps,dirichletsplit=False):
    ngmesh = ngm.Mesh()
    ngmesh.dim=3
    pnums = []

    for p in ngmeshbase.Points():
        x,y,z = p.p
        ngmesh.Add(ngm.MeshPoint(ngm.Pnt(x,y,0)))
    for i,p in enumerate(ngmeshbase.Points()):
        x,y,z = p.p
        pnums.append(ngmesh.Add(ngm.MeshPoint(ngm.Pnt(x,y,t_steps))))

    El1d = ngmeshbase.Elements1D()
    El2d = ngmeshbase.Elements2D()

    ngmesh.SetMaterial(1, "mat")
    for el in El2d:
        ngmesh.Add(ngm.Element3D(1, [
            pnums[el.points[0].nr-1]
            ,pnums[el.points[1].nr-1]
            ,pnums[el.points[2].nr-1]
            ,el.points[0].nr
            ,el.points[1].nr
            ,el.points[2].nr
            ]))

    fde = ngm.FaceDescriptor(surfnr=1,domin=1,bc=1)
    fde.bcname = "bottom"
    fdid = ngmesh.Add(fde)
    for el in El2d:
        ngmesh.Add(ngm.Element2D(fdid, [el.points[2].nr
            ,el.points[1].nr
            ,el.points[0].nr]))

    fde = ngm.FaceDescriptor(surfnr=2,domin=1,bc=2)
    fde.bcname = "top"
    fdid = ngmesh.Add(fde)
    for el in El2d:
        ngmesh.Add(ngm.Element2D(fdid, [pnums[el.points[0].nr-1]
            ,pnums[el.points[1].nr-1]
            ,pnums[el.points[2].nr-1]]))

    if dirichletsplit==False:
        fde = ngm.FaceDescriptor(surfnr=3,domin=1,bc=3)
        fde.bcname = "dirichlet"
        fdid = ngmesh.Add(fde)
        for el in El1d:
            ngmesh.Add(ngm.Element2D(fdid, [el.points[0].nr
                ,el.points[1].nr
                ,pnums[el.points[1].nr-1]
                ,pnums[el.points[0].nr-1]]))
            ngmesh.SetBCName(2,"dirichlet")

    if dirichletsplit==True:
        fde = ngm.FaceDescriptor(surfnr=3,domin=1,bc=3)
        fde.bcname = "front"
        fdid = ngmesh.Add(fde)
        for el in El1d:
            if ngmeshbase.Points()[el.points[0]][1]==0:
                ngmesh.Add(ngm.Element2D(fdid, [el.points[0].nr
                    ,el.points[1].nr
                    ,pnums[el.points[1].nr-1]
                    ,pnums[el.points[0].nr-1]]))

        fde = ngm.FaceDescriptor(surfnr=4,domin=1,bc=4)
        fde.bcname = "back"
        fdid = ngmesh.Add(fde)
        for el in El1d:
            if ngmeshbase.Points()[el.points[0]][1]==1:
                ngmesh.Add(ngm.Element2D(fdid, [el.points[0].nr
                    ,el.points[1].nr
                    ,pnums[el.points[1].nr-1]
                    ,pnums[el.points[0].nr-1]]))

        fde = ngm.FaceDescriptor(surfnr=5,domin=1,bc=5)
        fde.bcname = "left"
        fdid = ngmesh.Add(fde)
        for el in El1d:
            if ngmeshbase.Points()[el.points[0]][0]==0:
                ngmesh.Add(ngm.Element2D(fdid, [el.points[0].nr
                    ,el.points[1].nr
                    ,pnums[el.points[1].nr-1]
                    ,pnums[el.points[0].nr-1]]))

        fde = ngm.FaceDescriptor(surfnr=6,domin=1,bc=6)
        fde.bcname = "right"
        fdid = ngmesh.Add(fde)
        for el in El1d:
            if ngmeshbase.Points()[el.points[0]][0]==1:
                ngmesh.Add(ngm.Element2D(fdid, [el.points[0].nr
                    ,el.points[1].nr
                    ,pnums[el.points[1].nr-1]
                    ,pnums[el.points[0].nr-1]]))

        ngmesh.SetBCName(2,"front")
        ngmesh.SetBCName(3,"back")
        ngmesh.SetBCName(4,"left")
        ngmesh.SetBCName(5,"right")

    ngmesh.SetBCName(0,"bottom")
    ngmesh.SetBCName(1,"top")
    mesh = Mesh(ngmesh)
    return mesh

def CartSquare(N,t_steps,xscale=1,xshift=0,tscale=1):
    ngmesh = ngm.Mesh()
    ngmesh.SetGeometry(unit_square)
    ngmesh.dim = 2
    pnums = []
    for j in range(t_steps + 1):
        for i in range(N + 1):
            pnums.append(ngmesh.Add(ngm.MeshPoint(ngm.Pnt((i / N)*xscale-xshift, (j / t_steps)*tscale, 0))))

    foo = ngm.FaceDescriptor(surfnr=1,domin=1,bc=1)
    ngmesh.Add (foo)
    ngmesh.SetMaterial(1, "mat")
    for j in range(t_steps):
        for i in range(N):
            ngmesh.Add(ngm.Element2D(1, [pnums[i + j * (N + 1)],
                pnums[i + 1 + j * (N + 1)],
                pnums[i + 1 + (j + 1) * (N + 1)],
                pnums[i + (j + 1) * (N + 1)]]))

    fde = ngm.FaceDescriptor(surfnr=1,domin=1,bc=1)
    fde.bcname = "bottom"
    fdid = ngmesh.Add(fde)
    for i in range(N):
        ngmesh.Add(ngm.Element1D([pnums[i], pnums[i + 1]], index=1))

    fde = ngm.FaceDescriptor(surfnr=2,domin=1,bc=2)
    fde.bcname = "top"
    fdid = ngmesh.Add(fde)
    for i in range(N):
        ngmesh.Add(ngm.Element1D([pnums[i + t_steps * (N + 1)], pnums[i + 1 + t_steps * (N + 1)]], index=2))

    fde = ngm.FaceDescriptor(surfnr=3,domin=1,bc=3)
    fde.bcname = "dirichlet"
    fdid = ngmesh.Add(fde)
    for i in range(t_steps):
        ngmesh.Add(ngm.Element1D([pnums[N + i * (N + 1)], pnums[N + (i + 1) * (N + 1)]], index=3))
        ngmesh.Add(ngm.Element1D([pnums[0 + i * (N + 1)], pnums[0 + (i + 1) * (N + 1)]], index=3))


    ngmesh.SetBCName(0,"bottom")
    ngmesh.SetBCName(1,"top")
    ngmesh.SetBCName(2,"dirichlet")

    mesh = Mesh(ngmesh)
# print("boundaries" + str(mesh.GetBoundaries()))
    return mesh



def AddSurfElements2DBC(tpmesh,mesh1,mesh2):
    if mesh1.dim==2:
        ngm1 = mesh1.ngmesh;
        ngm2 = mesh2.ngmesh;
    else:
        ngm1 = mesh2.ngmesh;
        ngm2 = mesh1.ngmesh;
    els1 = ngm1.Elements2D()
    els2 = ngm2.Elements1D()
    # tpmesh.Add (FaceDescriptor(surfnr=1,domin=1,bc=1))

    fde = ngm.FaceDescriptor(surfnr=1,domin=1,bc=1)
    fde.bcname = "top"
    fdid = tpmesh.Add(fde)
    for elx in els1:
        vert_loc = elx.vertices
        vert_glob = []
        for vx in vert_loc:
            vert_glob.append(PointId((vx.nr-1)*len(ngm2.Points())+len(ngm2.Points())))
        tpmesh.Add(Element2D(fdid,vert_glob))

    fde = ngm.FaceDescriptor(surfnr=2,domin=1,bc=2)
    fde.bcname = "bottom"
    fdid = tpmesh.Add(fde)
    for elx in els1:
        vert_loc = elx.vertices
        vert_glob = []
        for vx in vert_loc:
            vert_glob.insert(0,PointId((vx.nr-1)*len(ngm2.Points())+1))
        tpmesh.Add(Element2D(fdid,vert_glob))

    els1 = ngm1.Elements1D()
    fde = ngm.FaceDescriptor(surfnr=3,domin=1,bc=3)
    fde.bcname = "dirichlet"
    fdid = tpmesh.Add(fde)
    for elx in els1:
        for ely in els2:
            vert_glob=[]
#            for vy in ely.vertices:
#                for vx in elx.vertices:
            vx = elx.vertices
            vy = ely.vertices
            vert_glob = [PointId((vx[1].nr-1)*len(ngm2.Points())+vy[0].nr),
                        PointId((vx[1].nr-1)*len(ngm2.Points())+vy[1].nr),
                        PointId((vx[0].nr-1)*len(ngm2.Points())+vy[1].nr),
                        PointId((vx[0].nr-1)*len(ngm2.Points())+vy[0].nr)]
            tpmesh.Add(Element2D(fdid,vert_glob))
    tpmesh.SetBCName(0,"top")
    tpmesh.SetBCName(1,"bottom")
    tpmesh.SetBCName(2,"dirichlet")
    return tpmesh


def TensorProdMesh(meshx,mesht):
    tpmesh = MakeMesh3D(meshx,mesht)
    return AddSurfElements2DBC(tpmesh,meshx,mesht)
