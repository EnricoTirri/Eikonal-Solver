/*
 * SimplexData.hpp
 *
 *  Created on: Jun 17, 2022
 *      Author: forma
 */

#ifndef EIKONAL_SIMPLEXDATA_HPP_
#define EIKONAL_SIMPLEXDATA_HPP_
#include "Eikonal_traits.hpp"
namespace Eikonal
{
template<std::size_t PHDIM>
  struct SimplexData
  {
	using AnisotropyM=typename Eikonal_traits<PHDIM>::AnisotropyM;
	using MMatrix=typename Eikonal_traits<PHDIM>::MMatrix;
	using Point=typename Eikonal_traits<PHDIM>::Point;

	//! This constructor just takes of vector that describes the simplex

	  SimplexData(std::array<std::array<double,PHDIM>,PHDIM+1u> const & p,
			  AnisotropyM const &M=AnisotropyM::Identity())
	  {
		  for (auto i=0u;i<PHDIM+1u;++i)
			  points[i]=Eigen::Map<Point>(const_cast<double*>(p[i].data()));
		  setup(M);
	  }

	  SimplexData(std::array<Point,PHDIM+1u> const & p,
			  AnisotropyM const &M=AnisotropyM::Identity()):
		  points{p}
	  {
		  setup(M);
	  };
	  std::array<Point,PHDIM+1u> points;
	  MMatrix MM_Matrix;
	  MMatrix E;
  private:
	  void setup(AnisotropyM const &M)
	  {
		  E.col(0)=points[PHDIM-1u]-points[0];//e13 or e11
		  if constexpr(PHDIM==3u)
    		{
			  E.col(1)=points[PHDIM-1u]-points[1]; //e23
    		}
		  E.col(PHDIM-1)=points[PHDIM]-points[PHDIM-1u]; //e34 or e23
		  MM_Matrix=E.transpose()*M*E;
	  }
  };
}




#endif /* EIKONAL_SIMPLEXDATA_HPP_ */
