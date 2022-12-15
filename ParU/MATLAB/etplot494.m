Parent = [ 8 68 68 68 67 67 67 67 67 67 67 67 67 67 67 67 66 66 66 66 66 66 66 28 28 28 27 28 66 66 66 66 65 65 65 65 65 65 65 65 64 64 64 64 64 64 64 64 64 63 63 63 63 63 63 63 62 62 62 62 62 62 63 64 65 66 67 68 -1];
Parent = Parent+1;
G = digraph(Parent(Parent~=0),find(Parent));
h = plot(G,'Layout','Layered', "ShowArrows", false);
h.Marker = 'o';
h.MarkerSize = 7;
