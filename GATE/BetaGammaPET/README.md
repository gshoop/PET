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
/gate/cylindricalPET/geometry/setRmin 16.0 cm
/gate/cylindricalPET/geometry/setRmax 20.2 cm
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
/gate/block/placement/setTranslation 18.0 0.0 0.0 cm
/gate/block/geometry/setXLength 4.0 cm
/gate/block/geometry/setYLength 4.0 cm
/gate/block/geometry/setZLength 15.0 cm
/gate/block/setMaterial Vacuum
/gate/block/vis/setVisible true
/gate/block/vis/setColor yellow
/gate/block/vis/forceWireFrame
```

# Geometry 2
![image](https://user-images.githubusercontent.com/44107373/219119650-5fcdc630-75a1-481e-8403-ac15a4994bf4.png)
![image](https://user-images.githubusercontent.com/44107373/219120103-4e2441df-bf20-42d6-b2c3-83f19981aa72.png)

Adding the 4x4x0.5 cm^3 CZT detector.
```
# PIXEL
/gate/block/daughters/name pixel
/gate/block/daughters/insert box
/gate/pixel/geometry/setXLength 4.0 cm
/gate/pixel/geometry/setYLength 4.0 cm
/gate/pixel/geometry/setZLength 0.5 cm
/gate/pixel/setMaterial CZT
/gate/pixel/setcolor green
/gate/pixel/forceSolid
/gate/pixel/setVisible true
```

# Geometry 3
![image](https://user-images.githubusercontent.com/44107373/219123212-c54fc40b-b814-46fc-9780-3262b0f8560a.png)

Repeating the CZT crystals in the z-direction for a total of 30.
```
# REPEATING PIXEL
/gate/pixel/repeaters/insert cubicArray
/gate/pixel/cubicArray/setRepeatNumberX 1
/gate/pixel/cubicArray/setRepeatNumberY 1
/gate/pixel/cubicArray/setRepeatNumberZ 30
/gate/pixel/cubicArray/setRepeatVector 0.0 0.0 0.5 cm
```

# Geometry 4
![image](https://user-images.githubusercontent.com/44107373/219128004-70959f9d-b0e9-40ad-bc85-4c7156a52cb0.png)

Repeating the CZT blocks in a ring for a total of 23 detector blocks around the circumference.
```
# RING REPEAT
/gate/block/repeaters/insert ring
/gate/block/ring/setRepeatNumber 23
```

```
# ATTACH VOLUMES AND SENSITIVE DETECTOR
/gate/systems/cylindricalPET/rsector/attach block
/gate/systems/cylindricalPET/module/attach pixel
/gate/pixel/attachCrystalSD
```

# Ring Diameter Optimization
The minimum ring size for the system that will use edge-on construction with 4x4x0.5 cm CZT crystals can be seen as a n-sided (n=23) regular polygon. The minimum ring size we can use can be calculated
from r = a/(2*tan(pi/n), where a=4 cm. This gives us a minimum ring radius of 14.55 cm. Many commercial brain/head PET systems have configurations with radii of 25 cm. In this case we made the ring radius
to be approximately 16 cm where we can see a few mm of spacing between each detector block in the ring.

# Sources

We would like to use Zirconium-89 (89Zr) for positron emission. 89Zr has a half-life of approximately 3.27 days or 78.48 hours. Zr decays by positron emission at a rate of 23% with average kinetic energy
of 395 keV. There are also gamma emissions at 99% with energy 909 keV.

