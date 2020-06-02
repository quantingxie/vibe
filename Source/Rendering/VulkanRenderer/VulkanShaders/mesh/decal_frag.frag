#version 460

layout (location = 0) out vec4 outputColor;
layout (location = 1) out vec4 outputSpecular;

layout (constant_id = 0) const uint numLights = 0;
layout (constant_id = 1) const bool tessellation = false;
layout (constant_id = 2) const bool shaded = true;
layout (constant_id = 3) const bool hasDiffuseTexture = false;
layout (constant_id = 4) const bool hasNormalTexture = false;
layout (constant_id = 5) const bool hasRoughnessTexture = false;
layout (constant_id = 6) const bool hasMetalnessTexture = false;
layout (constant_id = 7) const bool hasAmbientOcclusionTexture = false;
layout (constant_id = 8) const bool hasSubsurfaceScatteringTexture = false;
layout (constant_id = 9) const bool hasIrradianceCubemapTexture = false;
layout (constant_id = 10) const bool hasRadianceCubemapTexture = false;
layout (constant_id = 11) const bool hasBrdfLUTTexture = false;

struct light
{
    vec4 position;
    vec4 color;
    vec4 direction;
    ivec4 state; // 1 type, shadow map index
};

layout (set = 1, binding = 0) uniform globalUniforms
{
    mat4 inverseViewMatrices[2];
    mat4 inverseProjectionMatrices[2];
    vec4 resolution;
    light lights[16];
    mat4 lightMatrices[16];
} globals;

layout (set = 1, binding = 1) uniform localUniforms
{
    mat4 inverse[128];
    vec4 color;
    uint receivesShadows;
    float emissivity;
    float roughness;
    float metalness;
} locals;

layout (location = 0) in vertexData{
    vec3 position;
    vec3 normal;
    vec2 uv;
    mat3 TBN;
    vec3 cameraPosition;
    flat int index;
    flat uint view;
}vertex;

layout (set = 1, binding = 2) uniform sampler2D diffuseTexture;
layout (set = 1, binding = 3) uniform sampler2D normalTexture;
layout (set = 1, binding = 4) uniform sampler2D roughnessTexture;
layout (set = 1, binding = 5) uniform sampler2D metalnessTexture;
layout (set = 1, binding = 6) uniform sampler2D ambientOcclusionTexture;
layout (set = 1, binding = 7) uniform sampler2D subsurfaceScatteringTexture;
layout (set = 1, binding = 8) uniform sampler2DArray shadowArray;
layout (set = 1, binding = 9) uniform samplerCube irradianceCubemapTexture;
layout (set = 1, binding = 10) uniform samplerCube radianceCubemapTexture;
layout (set = 1, binding = 11) uniform sampler2D brdfLUTTexture;
layout (set = 1, binding = 12) uniform sampler2DArray aoBuffer;
layout (set = 1, binding = 13) uniform sampler2DArray depthAttachment;
layout (set = 1, binding = 14) uniform sampler2DArray normalAttachment;

// Constants
const float PI = 3.1415;

// Global variables
vec3 finalDiffuse = vec3(0);
vec3 finalSpecular = vec3(0);
vec3 diffuseColor = vec3(1, 1, 1);
vec3 normal = vec3(0, 0, 1);
float specularValue = 20;
vec3 cameraDirection = vec3(0, 0, 1);
float specularPow = 0;
vec3 specularColor = vec3(1, 1, 1);
float roughness = 1.0;
float metalness = 0.0;
float ambientOcclusion = 1.0;
float subsurfaceScattering = 0.0;
vec3 indirectDiffuse = vec3(0, 0, 0);
vec3 indirectSpecular = vec3(0, 0, 0);
float opacity = 1.0;
vec2 uv = vec2(0);
vec4 position = vec4(0, 0, 0, 1);

// Functions
void readTextures();
void calculateIndirectLighting();
void calculatePBRLighting(vec3 lightDirection, vec3 lightColor, float lightIntensity);
float geometryTerm(vec3 vector);
float squared(float x);
float calculateShadow(int index);
vec3 reconstructPosition(float depth);

void main(void)
{
    position = vec4(reconstructPosition(texture(depthAttachment, vec3(gl_FragCoord.xy/globals.resolution.xy, vertex.view)).x), 1);

    vec4 decalSpacePosition = locals.inverse[vertex.index] * position;
    uv = decalSpacePosition.xy + 0.5;

    if (abs(decalSpacePosition.x) > 0.5) discard;
    if (abs(decalSpacePosition.y) > 0.5) discard;
    if (abs(decalSpacePosition.z) > 0.5) discard;

    readTextures();

    normal = texture(normalAttachment, vec3(gl_FragCoord.xy/globals.resolution.xy, vertex.view)).rgb;

    // If it's 0, then there's a divide by zero error
    roughness = max(roughness * roughness, 0.0001);
    specularColor = mix(vec3(1.0), diffuseColor, metalness);

    cameraDirection = normalize(vertex.cameraPosition - position.xyz);

    vec3 lightRay = vec3(0);
    float lightIntensity = 0;

    for (int i = 0; i < numLights; i++)
    {
        int type = globals.lights[i].state.x;
        if (type > 1)
        {
            lightRay = position.xyz - globals.lights[i].position.xyz;
            float distanceSquared = dot(lightRay, lightRay);
            lightIntensity = globals.lights[i].color.a / distanceSquared;

            if (type == 3)
            {
                if (globals.lights[i].direction.a > dot(normalize(lightRay), normalize(globals.lights[i].direction.xyz)))
                {
                    lightIntensity = 0;
                }
            }
        }
        else
        {
            lightIntensity = globals.lights[i].color.a * calculateShadow(i);
            lightRay = globals.lights[i].direction.xyz;
        }

        calculatePBRLighting(normalize(lightRay), globals.lights[i].color.rgb, lightIntensity);
    }

    calculateIndirectLighting();

    finalDiffuse *= diffuseColor;
    outputColor = vec4(finalDiffuse + (diffuseColor * locals.emissivity), opacity);
    outputSpecular = vec4(finalSpecular, opacity);
}

void readTextures()
{
    diffuseColor = locals.color.rgb;
    if (hasDiffuseTexture)
    {
        vec4 tempColor = texture(diffuseTexture, uv).rgba;
        diffuseColor *= tempColor.rgb;
        opacity = tempColor.a;
    }

    if (hasNormalTexture)
    {
        normal = vertex.TBN * normalize((2.0 * texture(normalTexture, uv).rgb) - 1.0);
    }
    else
    {
        normal = vertex.TBN * normal;
    }

    if (hasRoughnessTexture)
    {
        roughness = texture(roughnessTexture, uv).r;
    }
    else
    {
        roughness = locals.roughness;
    }

    if (hasMetalnessTexture)
    {
        metalness = texture(metalnessTexture, uv).r;
    }
    else
    {
        metalness = locals.metalness;
    }

    if (hasSubsurfaceScatteringTexture)
    {
        subsurfaceScattering = texture(subsurfaceScatteringTexture, uv).r;
    }
}

void calculateIndirectLighting()
{
    // Fresnel term: Schlick's approximation
    vec3 F_0 = mix(vec3(0.04), diffuseColor, metalness);
    vec3 F = max((F_0) + (max(vec3(1.0 - roughness), F_0) - F_0) * pow(1.0 - max(dot(cameraDirection, normal), 0), 5), 0);

    vec3 reflection = reflect(-cameraDirection, normal);

    // Energy conservation
    vec3 k_s = F;
    vec3 k_d = (1 - k_s) * (1.0 - metalness);

    if (hasIrradianceCubemapTexture)
    {
        indirectDiffuse = texture(irradianceCubemapTexture, normal).rgb;
    }

    if (hasBrdfLUTTexture)
    {
        vec2 brdfValue = texture(brdfLUTTexture, vec2(max(dot(cameraDirection, normal), 0.0), roughness)).rg;

        if (hasRadianceCubemapTexture)
        {
            int numMipLevels = textureQueryLevels(radianceCubemapTexture);
            indirectSpecular = textureLod(radianceCubemapTexture, reflection, roughness * numMipLevels).rgb;
            indirectSpecular = indirectSpecular * (k_s * brdfValue.r + brdfValue.g);
        }
    }

    float ao = min(texture(aoBuffer, vec3(gl_FragCoord.xy / (textureSize(aoBuffer, 0) * 2).xy, vertex.view)).r, ambientOcclusion);
    finalDiffuse += indirectDiffuse * k_d * ao;
    finalSpecular += indirectSpecular * specularColor * ao;
}

void calculatePBRLighting(vec3 lightDirection, vec3 lightColor, float lightIntensity)
{
    // Lambert BRDF
    lightDirection *= -1;
    float diffusePow = 1.0 / PI;
    float l_dot_n = max(dot(normal, lightDirection), 0);

    vec3 halfway = normalize(normalize(lightDirection) + cameraDirection);

    // Cook-Torrance BRDF

    // Distribution term: Trowbridge-Reitz
    float roughness_squared = roughness * roughness;
    float D = roughness_squared / (PI * squared(squared(max(dot(halfway, normal), 0)) * (roughness_squared - 1) + 1.001));

    // Fresnel term: Schlick's approximation
    vec3 F_0 = mix(vec3(0.04), diffuseColor, metalness);
    vec3 F = max((F_0) + (1.0 - F_0) * pow(1.0 - max(dot(cameraDirection, halfway), 0), 5), 0);

    // Geometry term: Schlick's GGX
    float G = geometryTerm(cameraDirection) * geometryTerm(lightDirection);

    vec3 specularPow = (D * F * G) /
        (4 * max(dot(normal, cameraDirection), 0.001) * max(dot(normal, lightDirection), 0.001));

    // Energy conservation
    vec3 k_s = F;
    vec3 k_d = (1 - k_s) * (1.0 - metalness);

    finalDiffuse += k_d * diffusePow * lightColor * l_dot_n * lightIntensity;
    finalSpecular += specularPow * specularColor * lightColor * l_dot_n * lightIntensity;
}

float geometryTerm(vec3 vector)
{
    float k = roughness / 2;
    float vector_dot_n = max(dot(vector, normal), 0);
    return vector_dot_n / (vector_dot_n * (1 - k) + k);
}

float squared(float x)
{
    return x * x;
}

float calculateShadow(int index)
{
    int shadowMapIndex = globals.lights[index].state.y;
    if (locals.receivesShadows == 0 || shadowMapIndex == -1)
    {
        return 1;
    }

    float coverage = 0;
    vec4 pos = globals.lightMatrices[shadowMapIndex] * position;

    if (abs(pos.x) >= 1.0 || abs(pos.y) >= 1.0 || abs(pos.z) >= 1.0)
    {
        return 1;
    }

    pos.xy = 0.5 * pos.xy + 0.5;

    float resolution = globals.resolution.z;

    for (int x = -2; x < 3; x++)
    {
        for (int y = -2; y < 3; y++)
        {
            float depth = texture(shadowArray, vec3(pos.x + x / resolution, pos.y + y / resolution, shadowMapIndex)).r;

            if (pos.z > depth + 0.002)
            {                
                coverage += 1.0 / 25.0;
            }
        }
    }

    return (1 - coverage);
}

vec3 reconstructPosition(float depth)
{
    vec3 coords = vec3(gl_FragCoord.xy / globals.resolution.xy, 1);
    coords.x = 2 * coords.x - 1;
    coords.y = 2 * coords.y - 1;
    coords.z = depth;

    vec4 pos = vec4(coords, 1);
    pos = globals.inverseProjectionMatrices[vertex.view] * pos;
    pos /= pos.w;

    pos = globals.inverseViewMatrices[vertex.view] * pos;

    return pos.xyz;
}