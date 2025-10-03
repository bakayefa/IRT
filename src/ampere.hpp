#ifndef HYBRIDIR_AMPERE_HPP
#define HYBRIDIR_AMPERE_HPP

#include "vecfield.hpp"

#include <cstddef>
#include <iostream>

template<std::size_t dimension>
class Ampere
{
public:
    Ampere(std::shared_ptr<GridLayout<dimension>> grid)
        : m_grid{grid}
    {
        if (!m_grid)
            throw std::runtime_error("GridLayout is null");
    }

    void operator()(VecField<dimension> const& B, VecField<dimension>& J)
    {
        auto const dx = m_grid->cell_size(Direction::X);

        if constexpr (dimension == 1)
        {
            // TODO your code here
            auto const& Bx = B.x();
            auto const& By = B.y();
            auto const& Bz = B.z();
            
            auto& Jx = J.x();
            auto& Jy = J.y();
            auto& Jz = J.z();

            std::size_t N = Bx.data().size();  
            // For each grid point but Avoid first and last cells
            for (std::size_t ix = 1; ix < N-1; ++ix)            
            {
                Jx(ix) = 0.0;
                Jy(ix) = -(Bz(ix+1) - Bz(ix-1)) / (2*dx);
                Jz(ix) =  (By(ix+1) - By(ix-1)) / (2*dx);
            }

        }
        else
            throw std::runtime_error("Ampere not implemented for this dimension");
    }

private:
    std::shared_ptr<GridLayout<dimension>> m_grid;
};

#endif // HYBRIDIR_AMPERE_HPP
