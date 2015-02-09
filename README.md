# Physics-Simulator
Project From my college days. Initially built for the linux computers at UBC, I ported it to work on my windows machine.

This works on my machine with Visual Studios 2010 installed, and GLUT and openGL libraries linked.


This code has not been touched since I was in school (2013), so it may not be done amazingly well. But this was easily my favorite project from university, enough so that I added features to it after the class was over. The original was worked on with one of my classmates: Ignacio Rodriguez. The porting, bug fixing, and additional features was done by me.

I have all the Visual Studios project stuff up here, though you can rebuild it all yourself with just the main.cpp file.

I highly doubt I will ever touch this codebase again (even though it was tons of fun, and still needs a lot of work), but I was really proud of it. I could see myself taking this idea, and rebuilding it in some fashion in the future.

# What it does

This is a physics simulator. It can do electromagnetic and gravity forces, but primarily does gravity. Collisions can be turned on and off, and I put in some collision features like a backetball hoop, and a box. There are some preset scenarios (I think I gutted all the interesting ones, and the current ones are more debugging ones. sorry), and you can also "spawn" balls, with I believe the "m" key.

I found an old readme, which should explain most of what it does. This was written before the port/added features, but is probably the best documentation available.

# Old ReadMe

CJ Brassington
Ignacio Rodriguez


Our project is intended to create a realistic physics engine in openGl. It is based on particles, with given masses and charges. Our system is intended to realisitically model the forces that the particles exert on eachother due to charge and gravity. Also, our colision system is intended to realistically model the angle and force of the collision. I wont bother to explain the math behind them, as what I have figured out in the 15 or so times explaining it to other people, is that you either already know these things, or it takes a good 30 mins to explain.

To get the program to work, either click on "assmt" or in the command line type "./assmt"

When it opens, you will see a blank screen, this is because we decided to make a simulator and not a scene. Our physics simulator includes gravitational and coulomb forces, as well as dealing with collisions. The gravitaional force as well as the coulomb force is always on, however you can chose to turn collisions off. We based the code initialy off of A1, so we still have keyframes happening, though they are only there to make time go forward at this point.

We have 5 pre-made scenarios, and we also have several things that can be done to effect them. Pressing 'e' will make it earth like. What this means is that a green floor will appear and every object in the frame will be accelerated down at 9.8m/s. If/when the particles hit the ground, they will bounce back up. If you do not press 'e' then there will be no universal earth like gravity, so objects will only interact with each other. Pressing 'k' will turn collisions off. Collisions need to be on to have objects change their path when they hit into anything (the earth, another particle, a wall). Pressing 'W' (thats shift w) will turn on walls, which will form a box. Any particle hitting the walls of the box will be bounced back according to how they should. Pressing 'm' will draw a white particle of mass 1 and radius .25. This can be done as many times as you want, but be warned that when the animation starts, having more then 2 seems to usualy cause problems (see problems below). If you want to see where this particle will be placed before you do it, press 'M'. All these keys, with the exception of 'm', can be pressed again to reverse the effect.

When you get a scenario going, or to change where the 'm' balls will get placed, you have to move around using wasd as well as zx. wasd take care of forward/backward and left/right movement, and zx do up and down.

Scenario 1 was intended to look like a cool representation of a classical atom. It is in no way to scale, but the large nucleus is positively charged and the electrons are negatively charged. This is best viewed with collisions turned off.

Scenario 2 is a view of the solar system. This was initialy to scale, but you couldnt see anything but the sun and jupiter. Additionaly, even when we brought everything down by a factor of 10^11, the mass of the sun would not fit in a double, and making it any smaller would make the distances between the sun and jupiter much much less then 1. So we just approximated it, making the sun decently smaller, but keeping mercury, venus, earth and mars roughly proportional to each other. This should run well normaly, for some extra destructive fun put some "asteroids" (using m) in the path of a planet, and watch the consequences (might need a few hits as the planets all have masses in the thousands of kg's).

Scenario 3 was ment to show simple collisions. As is, it is fairly boring, but does show a collisions of 2 objects of the same mass. To spice it up, you can add particles for them to colide with, and put the walls up (move the camera in and up with w and x to be able to get a good view). Due to low speeds, yet high dt, it goes fairly fast if you put earth on.

Scenario 4 was ment to show off the earth feature. Really all it is is 3 bouncy balls, turn on 'e' and watch em go. As always, you can mess with this as you desire.

Scenario 5 was supposed to be the grand finale, but due to problems discusses below was a big let down. It was supposed to be around 20 balls that we would turn on earth and the walls, and watch the mayhem. But the balls would disapear whenever we go to the 5-10 particle mark. We left it in to show this in the demo. To see the problem for youself, get 5 going with 'e' on, it should fail very fast.



problems:
There are a few problems that we have long recognized, but cant seem to fix. The first is that while the math handles most collisions correctly, grazing collisions seem to unexpectedly add energy to the system (that is, the 2 objects both pick up speed). The math for collisions seems to be correct to me, and I have about 8 pages of scribbles supporting this, but none the less some small tweeks to it sometimes produce effects that solve the problem of grazing collisions, but then fail at direct ones. We decided to settle on what we have now, as it coveres the majority of what is in the simulator.

The second problem, which is a bigger one, and more of a mystery, is that the positions of particles will randomly becomes 'nan' in certain situations. We assumed we were either dividing by zero, or setting a double to something that wasnt a number somewhere in our code, but over the course of 2 hours, we at one point had, between if statements and setting of constants, everything done so that we could 100% guarentee it would not produce a 'nan' yet it still did. So we have given up, and just accept that this problem is there. It tends to pop up when you have too many objects, and seems to happen most often if 'e' is on. This really limited what we could show, as what we wanted scenario 5 to show was a bunch of bouncy balls going crazy as they interact with each other (in the box, and on earth of course).

KEYS:

'A' = start animation
spacebar = stop animation
1-5 = scenarios
wasd = moving the camera around
zx = moving camera up/down
e = Turn on/off a earth type setting
W = turn on/off the walls of a box
k = Turning off/on collisions
M = show a cursor showing where an added particle would go
m = add a particle
c = clear the simulator

