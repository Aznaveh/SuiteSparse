Parent = [11 2 11 9 8 8 7 8 9 10 11 12 13 -1 ];
Parent = Parent+1;
G = digraph(Parent(Parent~=0),find(Parent));
h = plot(G,'Layout','Layered', "ShowArrows", false);
h.Marker = 'o';
h.MarkerSize = 7;
