License
=========
Copyright (C) 2018 Jianxin Wang(jxwang@mail.csu.edu.cn),Chengqian Lu(chengqlu@csu.edu.cn)

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 3
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.

Jianxin Wang(jxwang@mail.csu.edu.cn),Chengqian Lu(chengqlu@csu.edu.cn)
School of Information Science and Engineering
Central South University
ChangSha
CHINA, 410083

Type: Package
Title: Predicting human lncRNA-disease associations based on geometric matrix completion
=================
Description: This package implements the GMCLDA algorithm with geometric matrix completion framework, predicting lncRNA-disease 
associations.

Files:
1.Dataset

1) lncRNA-disease.txt stores known lncRNA-disease association information;

2) lncRNA.txt and diseases.txt store lncRNA ids and disease ids, respectively;

2.Code
1) ADMM_PCG.m: function computing approximation using ADMM (Alternating Direction Method of Multipliers) 

2) gKernel.m: function computing Gaussian interaction profile kernel;

3) GMCLDA_demo.m: predict potential lncRNA-disease associations using geometric matrix completion; 

4) Laplacian.m : function computing Laplacian matrix;

5) prox_nuclearnorm.m: function calculating approximation of nuclearnorm;

6) sample_sparse_AtA.m: function sparse sampling

7) sample_sparse_t.m: function projecting observed elements

8) test_gamma.m: function teting the value of gamma

9) vec.m: function serializing 




