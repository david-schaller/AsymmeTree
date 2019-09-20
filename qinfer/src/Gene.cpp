#include "Gene.h"

namespace std
{
    template <>
    struct hash<Gene>
    {
        size_t operator()(const Gene& gene) const
        {
            // Return hash value of the identifier
            return hash<std::string>()(gene.getIdentifier());
        }
    };
}
