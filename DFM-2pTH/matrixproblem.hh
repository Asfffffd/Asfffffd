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
 * \brief The sub-problem for the matrix domain in the exercise on two-phase flow in fractured porous media.
 */
#ifndef DUMUX_COURSE_FRACTURESEXERCISE_MATRIX_PROBLEM_HH
#define DUMUX_COURSE_FRACTURESEXERCISE_MATRIX_PROBLEM_HH

// we need this in this test in order to define the domain
// id of the fracture problem (see function interiorBoundaryTypes())
#include <dune/common/indices.hh>

// include the base problem and properties we inherit from
#include <dumux/porousmediumflow/problem.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/boundarytypes.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/parallel/vectorcommdatahandle.hh>

#include <dumux/discretization/box/scvftoscvboundarytypes.hh>
#include <dumux/material/components/co2.hh>

#include "co2tables.hh"

namespace Dumux {

/*!
 * \ingroup MultiDomain
 * \ingroup MultiDomainFacet
 * \ingroup TwoPTests
 * \brief The sub-problem for the matrix domain in the exercise on two-phase flow in fractured porous media.
 */
template<class TypeTag>
class MatrixSubProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;

    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;
    using NumEqVector = GetPropType<TypeTag, Properties::NumEqVector>;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using PrimaryVariables = typename GridVariables::PrimaryVariables;
    using Scalar = typename GridVariables::Scalar;

    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using FluidState = GetPropType<TypeTag, Properties::FluidState>;

    using FVGridGeometry = typename GridVariables::GridGeometry;
    using FVElementGeometry = typename FVGridGeometry::LocalView;
    using SubControlVolume = typename FVGridGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVGridGeometry::SubControlVolumeFace;
    using GridView = typename FVGridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    using VolumeVariables = typename GridVariables::GridVolumeVariables::VolumeVariables;
    using ElementVolumeVariables = typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView;
    using ElementFluxVariablesCache = typename GridVariables::GridFluxVariablesCache::LocalView;



    // some indices for convenience
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    enum
    {
        pressureIdx = Indices::pressureIdx,
        saturationIdx = Indices::saturationIdx,
        temperatureIdx = Indices::temperatureIdx,

        //! Equation indices
		contiH2OEqIdx = Indices::conti0EqIdx,
        contiN2EqIdx = Indices::conti0EqIdx + FluidSystem::N2Idx,
        energyEqIdx = Indices::energyEqIdx,

        //! Phase indices
        wPhaseIdx = FluidSystem::H2OIdx,
        nPhaseIdx = FluidSystem::N2Idx,

        dimWorld = GridView::dimensionworld


    };

public:
    //! The constructor
    MatrixSubProblem(std::shared_ptr<const FVGridGeometry> fvGridGeometry,
                     std::shared_ptr<typename ParentType::SpatialParams> spatialParams,
                     const std::string& paramGroup = "")
    : ParentType(fvGridGeometry, spatialParams, paramGroup)
    , boundaryOverPressure_(getParamFromGroup<Scalar>(paramGroup, "Problem.BoundaryOverPressure"))
    , boundarySaturation_(getParamFromGroup<Scalar>(paramGroup, "Problem.BoundarySaturation"))
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
        values.setAllNeumann();

        if (globalPos[1] > this->gridGeometry().bBoxMax()[1] - 1e-6)
//        if (globalPos[1] < eps_)
            values.setAllDirichlet();

        return values;
    }

    //! Specifies the type of interior boundary condition to be used on a sub-control volume face
    BoundaryTypes interiorBoundaryTypes(const Element& element, const SubControlVolumeFace& scvf) const
    {
        BoundaryTypes values;

        // Here we set the type of condition to be used on faces that coincide
        // with a fracture. If Neumann is specified, a flux continuity condition
        // on the basis of the normal fracture permeability is evaluated. If this
        // permeability is lower than that of the matrix, this approach is able to
        // represent the resulting pressure jump across the fracture. If Dirichlet is set,
        // the pressure jump across the fracture is neglected and the pressure inside
        // the fracture is directly applied at the interface between fracture and matrix.
        // This assumption is justified for highly-permeable fractures, but lead to erroneous
        // results for low-permeable fractures.
        // Here, we consider "open" fractures for which we cannot define a normal permeability
        // and for which the pressure jump across the fracture is neglectable. Thus, we set
        // the interior boundary conditions to Dirichlet.
        // IMPORTANT: Note that you will never be asked to set any values at the interior boundaries!
        //            This simply chooses a different interface condition!
        // TODO dumux-course-task C
        // Change coupling conditions!
        values.setAllDirichlet();

        return values;
    }

    //! evaluates the Dirichlet boundary condition for a given position
    PrimaryVariables dirichletAtPos(const GlobalPosition& globalPos) const
    {
        // initialize values with the initial conditions
        auto values = initialAtPos(globalPos);

        // nitrogen is in contact with the domain on the center half of the lower boundary
        // TODO dumux-course-task A
        // Change boundary conditions and Dirichlet values!
//        if (globalPos[1] < 1e-6 && globalPos[0] > 25.0 && globalPos[0] < 75.0)
//            values[saturationIdx] = boundarySaturation_;
//        	values[temperatureIdx] = 283.0;

        return values;
    }

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

//    NumEqVector neumannAtPos (const GlobalPosition& globalPos) const
//    {
//    	NumEqVector values(0.0);
//    	FluidState fs;
//    	if(globalPos[0] < 75 && globalPos[0] > 25)
//		{
//            const auto initialValues = initialAtPos(globalPos);
//            fs.setPressure(nPhaseIdx, initialValues[pressureIdx]); // assume pressure equality here
//            fs.setPressure(wPhaseIdx, initialValues[pressureIdx]);
//            fs.setTemperature(nPhaseIdx,273.15);
//            fs.setTemperature(wPhaseIdx,273.15);
//
////			values[contiN2EqIdx] = injectionrate_ *FluidSystem::density(fs,nPhaseIdx);
//            values[contiN2EqIdx] = -1e-3; //kg/s/m2
////            values[energyEqIdx] = -darcyVelocity_*volVars.density(wPhaseIdx)
////                                     *IapwsH2O::liquidEnthalpy(temperatureHigh_, volVars.pressure(wPhaseIdx));
//			values[energyEqIdx] = values[contiN2EqIdx] * FluidSystem::heatCapacity(fs,nPhaseIdx) * 100 ;
//		}
//		return values;
//    }

    NumEqVector sourceAtPos(const GlobalPosition &globalPos) const
    {
        NumEqVector values(0.0);
        values = 0;
        return values;
    }

    NumEqVector neumannAtPos(const GlobalPosition &globalPos) const
    {
        NumEqVector values(0.0);

        if (globalPos[0] < 0.8 + eps_ && globalPos[0] > 0.2 - eps_)
        {
            // inject air. negative values mean injection
            values[contiN2EqIdx] = -1e-3; // kg/(s*m^2)

            // compute enthalpy flux associated with this injection [(J/(kg*s)]
            using FluidState = GetPropType<TypeTag, Properties::FluidState>;
            FluidState fs;

            const auto initialValues = initialAtPos(globalPos);
            fs.setPressure(wPhaseIdx, initialValues[pressureIdx]);
            fs.setPressure(nPhaseIdx, initialValues[pressureIdx]); // assume pressure equality here
            fs.setTemperature(wPhaseIdx, 273.15);
            fs.setTemperature(nPhaseIdx, 273.15);
//
//            // energy flux is mass flux times specific enthalpy
            values[energyEqIdx] = values[contiN2EqIdx]*FluidSystem::enthalpy(fs, nPhaseIdx);
        }

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

    Scalar boundaryOverPressure_;
    Scalar boundarySaturation_;
    static constexpr Scalar eps_ = 1.5e-7;
//    int injectionrate_;
};

} // end namespace Dumux

#endif
