Parent = [91 91 91 91 91 91 91 91 91 91 90 89 88 87 86 86 86 86 86 86 86 86 86 86 86 86 86 86 86 86 86 86 86 86 85 84 83 82 82 82 82 82 82 82 82 82 82 81 58 58 58 58 58 58 58 58 58 58 80 79 78 78 78 78 77 76 76 75 74 73 72 72 73 74 75 76 77 78 79 80 81 82 83 84 85 86 87 88 89 90 91 174 173 172 171 170 170 170 170 170 170 170 170 170 170 170 170 170 170 170 169 169 169 169 169 169 169 169 169 169 169 169 169 168 168 168 168 168 168 168 168 168 168 168 167 166 138 138 141 141 141 145 145 145 145 165 149 149 149 164 153 153 153 163 156 156 162 161 160 160 161 162 163 164 165 166 167 168 169 170 171 172 173 174 -1 
];
Parent = Parent+1;
G = digraph(Parent(Parent~=0),find(Parent));
h = plot(G,'Layout','Layered', "ShowArrows", false);
