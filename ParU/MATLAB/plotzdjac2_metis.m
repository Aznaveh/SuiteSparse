Parent = [4 3 3 4 5 6 7 8 9 10 11 12 -1 53 53 53 53 51 51 51 51 51 51 51 51 51 51 51 51 51 51 51 51 51 51 51 51 51 51 51 51 51 51 51 51 51 51 51 51 51 51 52 53 54 55 56 57 58 59 60 61 -1 63 64 65 66 67 68 69 135 129 128 128 128 128 128 128 128 128 80 128 127 127 127 127 127 127 127 127 127 127 127 127 127 127 126 126 126 126 126 126 126 126 126 126 126 126 126 126 126 126 126 126 126 126 126 126 125 125 125 125 125 125 125 125 126 127 128 129 130 131 132 133 134 135 136 137 138 139 140 141 142 270 178 178 178 149 149 149 178 177 177 177 176 176 176 176 176 176 176 176 176 176 176 176 176 176 176 176 176 176 176 176 176 176 176 177 178 179 180 181 182 183 184 185 186 187 188 224 191 191 215 214 214 213 213 213 213 212 212 212 212 212 212 212 212 212 212 212 212 212 212 213 214 215 216 217 218 219 220 221 222 223 224 269 228 228 228 232 231 231 232 246 240 240 240 240 239 239 240 241 242 243 244 245 246 268 248 250 250 257 256 256 254 255 256 257 258 259 260 261 262 263 264 265 266 267 268 269 270 570 570 273 277 277 277 277 302 297 297 297 297 297 297 296 296 296 296 296 296 296 296 296 296 296 296 297 298 299 300 301 302 303 304 305 349 343 342 342 342 342 342 342 342 342 342 342 342 342 342 342 342 342 342 342 342 342 342 342 342 342 342 342 342 342 342 342 342 342 342 342 342 343 344 345 346 347 348 349 438 352 352 355 355 355 361 358 358 360 360 361 397 389 389 389 389 389 389 388 388 388 388 388 388 388 388 388 388 388 388 388 388 388 388 388 388 388 388 389 390 391 392 393 394 395 396 397 437 399 432 428 428 428 428 428 427 426 426 426 426 426 426 426 426 426 426 426 426 426 426 426 426 426 426 426 426 427 428 429 430 431 432 433 434 435 436 437 438 439 569 481 480 480 480 480 478 478 478 478 478 478 478 478 478 478 478 478 478 459 478 478 478 478 478 478 478 478 478 478 478 478 478 478 478 478 478 478 478 479 480 481 482 483 484 485 486 487 488 519 492 491 492 496 495 495 496 497 498 499 500 518 503 503 513 511 511 511 510 509 510 511 512 513 514 515 516 517 518 519 568 522 522 529 526 525 526 527 528 529 530 531 566 536 535 535 536 565 539 539 551 543 543 543 544 545 546 547 548 549 550 551 564 557 556 556 556 557 562 560 560 561 562 563 564 565 566 567 568 569 570 -1 
];
Parent = Parent+1;
G = digraph(Parent(Parent~=0),find(Parent));
h = plot(G,'Layout','Layered', "ShowArrows", false);