// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
 /*!
  * \file
  * \ingroup MultiDomain
  * \ingroup MultiDomainFacet
  * \ingroup TwoPTests
  * \brief The sub-problem for the fracture domain in the exercise on two-phase flow in fractured porous media.
  */
#ifndef DUMUX_COURSE_FRACTURESEXERCISE_FRACTURE_PROBLEM_HH
#define DUMUX_COURSE_FRACTURESEXERCISE_FRACTURE_PROBLEM_HH

// include the base problem and properties we inherit from
#include <dumux/porousmediumflow/problem.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/boundarytypes.hh>

namespace Dumux {

/*!
 * \ingroup MultiDomain
 * \ingroup MultiDomainFacet
 * \ingroup TwoPTests
  * \brief The sub-problem for the fracture domain in the exercise on two-phase flow in fractured porous media.
 */
template<class TypeTag>
class FractureSubProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;

    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;
    using NumEqVector = GetPropType<TypeTag, Properties::NumEqVector>;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using ElementVolumeVariables = typename GridVariables::GridVolumeVariables::LocalView;
    using PrimaryVariables = typename GridVariables::PrimaryVariables;
    using Scalar = typename GridVariables::Scalar;

    using FVGridGeometry = typename GridVariables::GridGeometry;
    using FVElementGeometry = typename FVGridGeometry::LocalView;
    using SubControlVolume = typename FVGridGeometry::SubControlVolume;
    using GridView = typename FVGridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    // some indices for convenience
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    enum
    {
        pressureIdx = Indices::pressureIdx,
        saturationIdx = Indices::saturationIdx,
        temperatureIdx = Indices::temperatureIdx
    };

public:
    //! The constructor
    FractureSubProblem(std::shared_ptr<const FVGridGeometry> fvGridGeometry,
                       std::shared_ptr<typename ParentType::SpatialParams> spatialParams,
                       const std::string& paramGroup)
    : ParentType(fvGridGeometry, spatialParams, paramGroup)
    , aperture_(getParamFromGroup<Scalar>(paramGroup, "SpatialParams.Aperture"))
    {
        // initialize the fluid system, i.e. the tabulation
        // of water properties. Use the default p/T ranges.
        using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
        FluidSystem::init();
    }

    //! Specifies the type of boundary condition at a given position
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition& globalPos) const
    {
        BoundaryTypes values;

        // We only use no-flow boundary conditions for all immersed fractures
        // in the domain (fracture tips that do not touch the domain boundary)
        // Otherwise, we would lose mass leaving across the fracture tips.
        values.setAllNeumann();

        // However, there is one fracture reaching the top boundary. For this
        // fracture tip we set Dirichlet Bcs as in the matrix domain
        // TODO dumux-course-task A
        // Change boundary conditions
        if (globalPos[1] > this->gridGeometry().bBoxMax()[1] - 1e-6)
            values.setAllDirichlet();

        return values;
    }

    //! Evaluate the source term in a sub-control volume of an element
    NumEqVector source(const Element& element,
                       const FVElementGeometry& fvGeometry,
                       const ElementVolumeVariables& elemVolVars,
                       const SubControlVolume& scv) const
    {
        // evaluate sources from bulk domain using the function in the coupling manager
        auto source = couplingManagerPtr_->evalSourcesFromBulk(element, fvGeometry, elemVolVars, scv);

        // these sources are in kg/s, divide by volume and extrusion to have it in kg/s/m³
        source /= scv.volume()*elemVolVars[scv].extrusionFactor();
        return source;
    }

    //! Set the aperture as extrusion factor.
    Scalar extrusionFactorAtPos(const GlobalPosition& globalPos) const
    {
        // We treat the fractures as lower-dimensional in the grid,
        // but we have to give it the aperture as extrusion factor
        // such that the dimensions are correct in the end.
        return aperture_;
    }

    //! evaluates the Dirichlet boundary condition for a given position
    PrimaryVariables dirichletAtPos(const GlobalPosition& globalPos) const
    { return initialAtPos(globalPos); }

    //! evaluate the initial conditions
    PrimaryVariables initialAtPos(const GlobalPosition& globalPos) const
    {
        // For the grid used here, the height of the domain is equal
        // to the maximum y-coordinate
        const auto domainHeight = this->gridGeometry().bBoxMax()[1];

        // we assume a constant water density of 1000 for initial conditions!
        const auto& g = this->spatialParams().gravity(globalPos);
        PrimaryVariables values;
        Scalar densityW = 1000.0;
        values[pressureIdx] = 1e5 - (domainHeight + 6000 - globalPos[1])*densityW*g[1];
        values[saturationIdx] = 0.0;
        values[temperatureIdx] = 383.0 + (domainHeight - globalPos[1])*0.03;

//        if (globalPos[0] > 20 - eps_ && globalPos[0] < 50 + eps_ && globalPos[1] > 5 - eps_ && globalPos[1] < 50 + eps_)
//            values[temperatureIdx] = 283;

        return values;
    }

    //! returns the temperature in \f$\mathrm{[K]}\f$ in the domain
//    Scalar temperature() const
//    { return 283.15; /*10°*/ }

    //! sets the pointer to the coupling manager.
    void setCouplingManager(std::shared_ptr<CouplingManager> cm)
    { couplingManagerPtr_ = cm; }

    //! returns reference to the coupling manager.
    const CouplingManager& couplingManager() const
    { return *couplingManagerPtr_; }

private:
    std::shared_ptr<CouplingManager> couplingManagerPtr_;
    Scalar aperture_;
    bool isExercisePartA_;
    bool isExercisePartB_;
    bool isExercisePartC_;
    static constexpr Scalar eps_ = 1.5e-7;
};

} // end namespace Dumux

#endif
