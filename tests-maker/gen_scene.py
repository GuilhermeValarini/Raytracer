import random
from datetime import datetime

random.seed(datetime.now())

sceneResName = ["0_hd", "1_full_hd", "2_4k"]
sceneRes = [[1280, 720], [1920, 1080], [3840, 2160]]
positionBorder = [[-10, 10], [-5, 5], [-100, -20]]
radiusBorder = [0, 3]
colorBorder= [[0, 1], [0, 1], [0, 1]]
reflectionBorder = [0, 1]

def createScene(spheres, scene):
    spheresText = "";
    for sphere in range(spheres):
        pos = []
        col = []
        for coord in range(3):
            pos.append(random.uniform(positionBorder[coord][0], positionBorder[coord][1]))
            col.append(random.uniform(colorBorder[coord][0], colorBorder[coord][1]))
        rad = random.uniform(radiusBorder[0], radiusBorder[1])
        ref = random.uniform(reflectionBorder[0], reflectionBorder[1])
        spheresText += "{:.2f} {:.2f} {:.2f} {:.2f} {:.2f} {:.2f} {:.2f} {:.2f}\n".format(pos[0], pos[1], pos[2], rad, col[0], col[1], col[2], ref)
    spheresText += "{:.2f} {:.2f} {:.2f} {:.2f} {:.2f} {:.2f} {:.2f} {:.2f}\n".format(0, 20, -30, 3, 0, 0, 0, 0, 3)

    for res, resName in zip(sceneRes, sceneResName):
        fileText = "{} {}\n".format(res[0], res[1])
        fileText += "{} 1\n".format(spheres)
        fileText += spheresText
        with open("scene{}_{}.in".format(scene, resName), "w") as testScene:
            testScene.write(fileText);

def createRandom():
    numberScenes = 5
    maxSpheres = 50
    for scene in range(numberScenes):
        spheres = random.randrange(1, maxSpheres)
        createScene(spheres, scene)

def create(spheres):
    for scene, sphere in zip(range(len(spheres)), spheres):
        createScene(sphere, scene)

numSpheres = [10, 40, 70, 100]
create(numSpheres)
