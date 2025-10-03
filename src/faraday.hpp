#ifndef HYBRIDIR_FARADAY_HPP
#define HYBRIDIR_FARADAY_HPP

#include "vecfield.hpp"
#include "gridlayout.hpp"
#include "utils.hpp"

#include <cstddef>
#include <iostream>
#include <memory>

template<std::size_t dimension>
class Faraday
{
    // TODO implement the Faraday class, hint - get inspiration from Ampere
public:
    Faraday(std::shared_ptr<GridLayout<dimension>> grid)
        : m_grid{grid}
    {
        if (!m_grid)
            throw std::runtime_error("GridLayout is null");
    }

    void operator()(VecField<dimension> const& E, VecField<dimension> const& B, VecField<dimension> Bnew)
    {
        auto const dx = m_grid->cell_size(Direction::X);

        if constexpr (dimension == 1)
        {
            // TODO your code here
            auto const& Ex = E.x();
            auto const& Ey = E.y();
            auto const& Ez = E.z();
            
            auto const& Bx = B.x();
            auto const& By = B.y();
            auto const& Bz = B.z();

            auto& Bnew_x = Bnew.x();
            auto& Bnew_y = Bnew.y();
            auto& Bnew_z = Bnew.z();

            std::size_t N = Bx.data().size();  
            // For each grid point but Avoid first and last cells
            for (std::size_t ix = 1; ix < N-1; ++ix)            
            {
                Bnew_x(ix) = Bx(ix);
                Bnew_y(ix) = By(ix) + (Ez(ix+1) - Ez(ix-1))/ (2.0 * dx);
                Bnew_z(ix) = Bz(ix) - (Ey(ix+1) - Ey(ix-1))/ (2.0 * dx);
            }

        }
        else
            throw std::runtime_error("Faraday not implemented for this dimension");
    }
private:
    std::shared_ptr<GridLayout<dimension>> m_grid;
};


#endif // HYBRIDIR_FARADAY_HPP
