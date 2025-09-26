#ifndef PHARE_CORE_NUMERICS_PUSHER_PUSHER_HPP
#define PHARE_CORE_NUMERICS_PUSHER_PUSHER_HPP


#include "vecfield.hpp"
#include "particle.hpp"

#include <cstddef>
#include <vector>


template<std::size_t dimension>
class Pusher //every Pusher knows about the grid it’s working on and the size of the time step it advances particles by.
{
protected:
    std::shared_ptr<GridLayout<dimension>> layout_;
    double dt_;

public:
    Pusher(std::shared_ptr<GridLayout<dimension>> layout, double dt)
        : layout_(layout)
        , dt_(dt)
    {
    }

    virtual void operator()(std::vector<Particle<dimension>>& particles,
                            VecField<dimension> const& E, VecField<dimension> const& B)
        = 0;

    virtual ~Pusher() {}
};

// Pusher<3> myPusher(layout, 0.01); would make a Pusher in 3D with timestep 0.01.


template<std::size_t dimension>
class Boris : public Pusher<dimension>
{
public:
    Boris(std::shared_ptr<GridLayout<dimension>> layout, double dt)
        : Pusher<dimension>{layout, dt}
    {
    }

    void operator()(std::vector<Particle<dimension>>& particles, VecField<dimension> const& E,
                    VecField<dimension> const& B) override
    {
        for (auto& particle : particles)
        {
            // TODO implement the Boris pusher

            auto dx = this->layout_->cell_size(Direction::X); 
            auto dt = this->dt_;

            
            particle.position[0] += particle.v[0]*dt;
                    
            
            auto iCell = static_cast<int>(particle.position[0] / dx); // (taken from layout mesh size) 
            auto remainder = (particle.position[0] / dx) - iCell;

            auto Ex = interpolate(E.x, iCell, remainder);
            auto Ey = interpolate(E.y, iCell, remainder);
            auto Ez = interpolate(E.z, iCell, remainder);

            auto Bx = interpolate(B.x, iCell, remainder);
            auto By = interpolate(B.y, iCell, remainder);
            auto Bz = interpolate(B.z, iCell, remainder);

            auto factor = (particle.charge*dt)/(2*particle.mass);
            // auto vminus = particle.velocity[0] + factor*E.x;
            // auto t = ((particle.charge*dt)/(2*particle.mass))*B;

            auto vminus = particle.v[0] + factor * Ex;
            auto vplus  = vminus;  
            
            particle.v[0] = vplus + factor*Ex;
            particle.position[0] += particle.v[0] * (dt/2);

        }
    }

private:
    double interpolate(Field<dimension> const& field, int iCell, double reminder) const
    {
        if (this->layout_->centerings(field.quantity())[0] == this->layout_->dual)
        {
            if (reminder < 0.5)
                return field(iCell - 1) * (1.0 - reminder) + field(iCell) * reminder;
            else
                return field(iCell) * (1.0 - reminder) + field(iCell + 1) * reminder;
        }
        return field(iCell) * (1.0 - reminder) + field(iCell + 1) * reminder;
    }
};


#endif
