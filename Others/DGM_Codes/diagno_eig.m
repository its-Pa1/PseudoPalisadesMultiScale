mat= [0.358534,0.00191345;0.00191345,0.306411];
[eigvec,eigval]=eig(mat);
transpose(eigvec)*mat*eigvec;
