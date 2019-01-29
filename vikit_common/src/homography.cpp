/*
 * homography.cpp
 * Adaptation of PTAM-GPL HomographyInit class.
 * https://github.com/Oxford-PTAM/PTAM-GPL
 * Licence: GPLv3
 * Copyright 2008 Isis Innovation Limited
 *
 *  Created on: Sep 2, 2012
 *      by: cforster
 */

#include <vikit/homography.h>
#include <opencv2/opencv.hpp>
#include <stdio.h>

namespace vk {

Homography::
Homography(const vector<Vector2d, aligned_allocator<Vector2d> >& _fts1,
           const vector<Vector2d, aligned_allocator<Vector2d> >& _fts2,
           double _error_multiplier2,
           double _thresh_in_px) :
   thresh(_thresh_in_px),
   error_multiplier2(_error_multiplier2),
   fts_c1(_fts1),
   fts_c2(_fts2)
{
}

void Homography::
calcFromPlaneParams(const Vector3d& n_c1, const Vector3d& xyz_c1)
{
  double d = n_c1.dot(xyz_c1); // normal distance from plane to KF
  H_c2_from_c1 = T_c2_from_c1.rotation_matrix() + (T_c2_from_c1.translation()*n_c1.transpose())/d;
}

/**
 * @brief 通过两幅图像中的匹配关系,调用OpenCV函数来得到单应矩阵
 * 
 */
void Homography::
calcFromMatches()
{
  vector<cv::Point2f> src_pts(fts_c1.size()), dst_pts(fts_c1.size());
  for(size_t i=0; i<fts_c1.size(); ++i)
  {
    src_pts[i] = cv::Point2f(fts_c1[i][0], fts_c1[i][1]);
    dst_pts[i] = cv::Point2f(fts_c2[i][0], fts_c2[i][1]);
  }

  // TODO: replace this function to remove dependency from opencv!
  cv::Mat cvH = cv::findHomography(src_pts, dst_pts, CV_RANSAC, 2./error_multiplier2);
  H_c2_from_c1(0,0) = cvH.at<double>(0,0);
  H_c2_from_c1(0,1) = cvH.at<double>(0,1);
  H_c2_from_c1(0,2) = cvH.at<double>(0,2);
  H_c2_from_c1(1,0) = cvH.at<double>(1,0);
  H_c2_from_c1(1,1) = cvH.at<double>(1,1);
  H_c2_from_c1(1,2) = cvH.at<double>(1,2);
  H_c2_from_c1(2,0) = cvH.at<double>(2,0);
  H_c2_from_c1(2,1) = cvH.at<double>(2,1);
  H_c2_from_c1(2,2) = cvH.at<double>(2,2);
}

size_t Homography::
computeMatchesInliers()
{
  inliers.clear(); inliers.resize(fts_c1.size());
  size_t n_inliers = 0;
  for(size_t i=0; i<fts_c1.size(); i++)
  {
    Vector2d projected = project2d(H_c2_from_c1 * unproject2d(fts_c1[i]));
    Vector2d e = fts_c2[i] - projected;
    double e_px = error_multiplier2 * e.norm();
    inliers[i] = (e_px < thresh);
    n_inliers += inliers[i];
  }
  return n_inliers;

}

bool Homography::
computeSE3fromMatches()
{
  calcFromMatches();
  bool res = decompose();
  if(!res)
    return false;
  computeMatchesInliers();
  findBestDecomposition();
  T_c2_from_c1 = decompositions.front().T;
  return true;
}


//HERE

/**
 * @brief 从单应矩阵恢复相机位姿RT
 * 
 * @return true   恢复成功
 * @return false  恢复失败
 */
bool Homography::
decompose()
{
  
  //不知道清空了什么 TODO
  decompositions.clear();
  //对计算出来的单应矩阵 H_c2_from_c1 进行求解,这里是直接对这个矩阵进行了奇异值分解,并且计算出来了两个奇异向量矩阵 UV
  JacobiSVD<MatrixXd> svd(H_c2_from_c1, ComputeThinU | ComputeThinV);
  //得到奇异值
  Vector3d singular_values = svd.singularValues();
  //提取奇异值
  double d1 = fabs(singular_values[0]); // The paper suggests the square of these (e.g. the evalues of AAT)
  double d2 = fabs(singular_values[1]); // should be used, but this is wrong. c.f. Faugeras' book.
  double d3 = fabs(singular_values[2]); // TODO 上面的两行注释中是什么意思?
  //得到两个奇异向量矩阵
  Matrix3d U = svd.matrixU();
  Matrix3d V = svd.matrixV();                    // VT^T
  //后面就是算法的正常工作了
  double s = U.determinant() * V.determinant();
  //注意这里的后缀,P=postive,M=???不知道,但是在Latex中符号\pm是表示正负号的意思,所以这里表示存储d2的数值没有考虑正负号
  double dPrime_PM = d2;

  //NOTICE  和ORB不同,SVO中求解的话,不考虑d是否正负,但是只会在d1!=d2!=d3的情况下才会进行求解
  int nCase;
  if(d1 != d2 && d2 != d3)
    nCase = 1;
  else if( d1 == d2 && d2 == d3)
    nCase = 3;
  else
    nCase = 2;

  if(nCase != 1)
  {
    printf("FATAL Homography Initialization: This motion case is not implemented or is degenerate. Try again. ");
    return false;
  }

  //同样的,x1 x3前面有个正负符号,这里也是用于表示这里存储的数值还需要经过一个正负符号的处理
  double x1_PM;
  double x2;
  double x3_PM;

  // All below deals with the case = 1 case.
  // Case 1 implies (d1 != d3)
  { // Eq. 12
  //吴博师兄ppt中的式17,暂时没有考虑正负号
    x1_PM = sqrt((d1*d1 - d2*d2) / (d1*d1 - d3*d3));
    x2    = 0;
    x3_PM = sqrt((d2*d2 - d3*d3) / (d1*d1 - d3*d3));
  };
  //然后把正负号在这里考虑
  double e1[4] = {1.0,-1.0, 1.0,-1.0};
  double e3[4] = {1.0, 1.0,-1.0,-1.0};

  Vector3d np;
  //这个类?结构体? 中存储了分解H矩阵所得到的旋转矩阵,平移向量和世界平面的法向量
  HomographyDecomposition decomp;

  //然后现在才开始就d是否大于0展开分类讨论
  // Case 1, d' > 0:
  decomp.d = s * dPrime_PM;
  //对于所有可能的符号组合展开遍历
  for(size_t signs=0; signs<4; signs++)
  {
    // Eq 13
    //吴博师兄PPT 式 18 19 , 计算旋转矩阵R'.
    //NOTE 这里只是暂时存储, 在这个函数的最后面还要进行最后一步解算,从R'->R
    decomp.R = Matrix3d::Identity();
    //在这里就能够看到对x1_PM 和 x3_PM 所进行的正负号处理了
    double dSinTheta = (d1 - d3) * x1_PM * x3_PM * e1[signs] * e3[signs] / d2;
    double dCosTheta = (d1 * x3_PM * x3_PM + d3 * x1_PM * x1_PM) / d2;
    decomp.R(0,0) = dCosTheta;
    decomp.R(0,2) = -dSinTheta;
    decomp.R(2,0) = dSinTheta;
    decomp.R(2,2) = dCosTheta;

    // Eq 14
    //计算平移向量t', 吴博师兄PPT  式 20
    decomp.t[0] = (d1 - d3) * x1_PM * e1[signs];
    decomp.t[1] = 0.0;
    decomp.t[2] = (d1 - d3) * -x3_PM * e3[signs];

    //吴博师兄PPT 式17 计算世界平面的法向量
    np[0] = x1_PM * e1[signs];
    np[1] = x2;
    np[2] = x3_PM * e3[signs];
    //对于它倒是直接进行解算了
    decomp.n = V * np;

    decompositions.push_back(decomp);
  }

  // Case 1, d' < 0:
  //对于d<0的操作也是相似的,为了省时间下面就暂时跳过去了
  decomp.d = s * -dPrime_PM;
  for(size_t signs=0; signs<4; signs++)
  {
    // Eq 15
    decomp.R = -1 * Matrix3d::Identity();
    double dSinPhi = (d1 + d3) * x1_PM * x3_PM * e1[signs] * e3[signs] / d2;
    double dCosPhi = (d3 * x1_PM * x1_PM - d1 * x3_PM * x3_PM) / d2;
    decomp.R(0,0) = dCosPhi;
    decomp.R(0,2) = dSinPhi;
    decomp.R(2,0) = dSinPhi;
    decomp.R(2,2) = -dCosPhi;

    // Eq 16
    decomp.t[0] = (d1 + d3) * x1_PM * e1[signs];
    decomp.t[1] = 0.0;
    decomp.t[2] = (d1 + d3) * x3_PM * e3[signs];

    np[0] = x1_PM * e1[signs];
    np[1] = x2;
    np[2] = x3_PM * e3[signs];
    decomp.n = V * np;

    decompositions.push_back(decomp);
  }

  // Save rotation and translation of the decomposition
  // 对前面得到的R' t'进行解算,直接得到 R t
  for(unsigned int i=0; i<decompositions.size(); i++)
  {
    Matrix3d R = s * U * decompositions[i].R * V.transpose();
    //TODO 但是这里可以说是没有进行过任何的处理哈,直接就使用了解算之后的t了,也没有进行任何的尺度固定处理
    //猜测可能是和ORB和SVO的目的不一样 ,SVO主要就是使用单应矩阵来恢复相机的运动,而只有第一帧进行初始化的时候
    //才需要进行单位化
    Vector3d t = U * decompositions[i].t;
    decompositions[i].T = Sophus::SE3(R, t);
  }
  return true;
}

bool operator<(const HomographyDecomposition lhs, const HomographyDecomposition rhs)
{
  return lhs.score < rhs.score;
}

/**
 * @brief 寻找到最佳的单应矩阵
 * 
 */
void Homography::
findBestDecomposition()
{
  //一般地, 这里应该有八和单应矩阵需要进行分解计算
  assert(decompositions.size() == 8);

  //step 1 深度同号约束
  for(size_t i=0; i<decompositions.size(); i++)
  {
    HomographyDecomposition &decom = decompositions[i];
    size_t nPositive = 0;
    //遍历参考帧中的所有的特征点
    for(size_t m=0; m<fts_c1.size(); m++)
    {
      //只对内点进行操作
      if(!inliers[m])
        continue;
      const Vector2d& v2 = fts_c1[m];
      //计算x2=H*x1的深度，并且判断空间点在两个帧中的深度是否同号
      double dVisibilityTest = (H_c2_from_c1(2,0) * v2[0] + H_c2_from_c1(2,1) * v2[1] + H_c2_from_c1(2,2)) / decom.d;
      if(dVisibilityTest > 0.0)
        nPositive++;
    }
    //排序用,符号是因为下面的排序函数是从小到大
    decom.score = -nPositive;
  }

  //根据刚才满足"深度同号"约束的点的个数,对八个单应矩阵进行排序
  sort(decompositions.begin(), decompositions.end());
  //由于根据计算结果,其中有四个的单应矩阵计算的得到的点深度比几乎都为负数,所以这里调整需要求解的总长度为四
  decompositions.resize(4);

  //step 2 平面深度为正
  //遍历剩下的四个单应矩阵
  for(size_t i=0; i<decompositions.size(); i++)
  {
    HomographyDecomposition &decom = decompositions[i];
    int nPositive = 0;
    //对于参考帧中的每一个点
    for(size_t m=0; m<fts_c1.size(); m++)
    {
      if(!inliers[m])
        continue;
      //将这个特征点反投影到空间中
      Vector3d v3 = unproject2d(fts_c1[m]);
      // \f$ \mathbf{n}^{\text{T}} \cdot \mathbf{X} / d = 1  \f$
      double dVisibilityTest = v3.dot(decom.n) / decom.d;
      //由于外点的存在,使得即使是对于正确的单应矩阵,也不能够保证所有的点都可以得到1,所以这里就只统计大于0的数目
      //TODO 但是为什么可以这样操作呢?
      if(dVisibilityTest > 0.0)
        nPositive++;
    };
    decom.score = -nPositive;
  }

  //同样选取深度最能够满足要求的点,选最好的前两个
  sort(decompositions.begin(), decompositions.end());
  decompositions.resize(2);

  // According to Faugeras and Lustman, ambiguity exists if the two scores are equal
  // but in practive, better to look at the ratio!
  //step 3 在这里作者觉得直接计算比例的话比较靠谱,如果这个比例超过0.9那么就直接选大的
  double dRatio = (double) decompositions[1].score / (double) decompositions[0].score;

  if(dRatio < 0.9) // no ambiguity!
    decompositions.erase(decompositions.begin() + 1);
  else  // two-way ambiguity. Resolve by sampsonus score of all points.
  {
    //误差的极限值? TODO 
    double dErrorSquaredLimit  = thresh * thresh * 4;
    //分别计算两个sampsonus的评分
    double adSampsonusScores[2];
    for(size_t i=0; i<2; i++)
    {
      Sophus::SE3 T = decompositions[i].T;
      //计算得到本征矩阵
      Matrix3d Essential = T.rotation_matrix() * sqew(T.translation());
      //累计在这个本征矩阵下,所有特征点到重投影极线的距离
      double dSumError = 0;
      for(size_t m=0; m < fts_c1.size(); m++ )
      {
        //对当前遍历到的点可计算sapsonus误差
        double d = sampsonusError(fts_c1[m],    //参考帧中的特征点
                                  Essential,    //上面的得到的本征矩阵
                                  fts_c2[m]);   //当前帧中的特征点
        //避免极端外点的影响
        if(d > dErrorSquaredLimit)
          d = dErrorSquaredLimit;
        dSumError += d;
      }
      adSampsonusScores[i] = dSumError;
    }

    if(adSampsonusScores[0] <= adSampsonusScores[1])
      decompositions.erase(decompositions.begin() + 1);
    else
      decompositions.erase(decompositions.begin());
  }
}


} /* end namespace vk */

