# Geometry 1
![image](https://user-images.githubusercontent.com/44107373/219108093-d059239e-7999-4801-b793-4c7d6460dd13.png)

This geometry was created with the following.

```
# WORLD GEOMETRY
/gate/world/geometry/setXLength 60. cm
/gate/world/geometry/setYLength 60. cm
/gate/world/geometry/setZLength 40. cm
/gate/world/setMaterial Vacuum
/gate/world/vis/setVisible true
/gate/world/vis/setColor white
/gate/world/vis/forceWireFrame
```

```
# CYLINDRICAL PET
/gate/world/daughters/name cylindricalPET
/gate/world/daughters/insert cylinder
/gate/cylindricalPET/placement/setTranslation 0.0 0.0 0.0 cm
/gate/cylindricalPET/geometry/setRmin 20.0 cm
/gate/cylindricalPET/geometry/setRmax 25.0 cm
/gate/cylindricalPET/geometry/setHeight 15.0 cm
/gate/cylindricalPET/setMaterial Vacuum
/gate/cylindricalPET/vis/setVisible true
/gate/cylindricalPET/vis/setColor cyan
/gate/cylindricalPET/vis/forceWireFrame
```

```
# DETECTOR BLOCK
/gate/cylindricalPET/daughters/name block
/gate/cylindricalPET/daughters/insert box
/gate/block/placement/setTranslation 22.2 0.0 0.0 cm
/gate/block/geometry/setXLength 4.0 cm
/gate/block/geometry/setYLength 4.0 cm
/gate/block/geometry/setZLength 15.0 cm
/gate/block/setMaterial Vacuum
/gate/block/vis/setVisible true
/gate/block/vis/setColor yellow
/gate/block/vis/forceWireFrame
```
