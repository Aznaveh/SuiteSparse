Parent = [2 2 34 33 33 9 8 8 9 10 32 30 29 18 17 16 17 18 28 20 28 27 27 27 25 26 27 28 29 30 31 32 33 34 35 -1 ];
Parent = Parent+1;
G = digraph(Parent(Parent~=0),find(Parent));
h = plot(G,'Layout','Layered', "ShowArrows", false);
