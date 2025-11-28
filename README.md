# Project Announcement

This project based on the following four components:

*   **Mid-term written report** (Due: Dec 3, team-based, up to 5 points)
*   **Mid-term oral presentation** (Temporarily scheduled on Dec 4, team-based, up to 5 points)
*   **In-class presentation** (Dec 31, team-based, up to 15 points)
*   **Final written report** (Due: Jan 14, **individual**-based, up to 15 points, hard deadline, no extensions)

The total project score is 40+bonus points.

## 1 Mid-term Written Report

Submit a report summarizing your project progress using the SIGGRAPH template. One **full** page is sufficient for this report, but do not be more than two pages. You need to include the following information:

*   The project topic you chose
*   Your project goal and the technical points you intended to include
*   A detailed schedule for the whole project
*   The technical aspects you have finished, including methods used and any visual results (if available)
*   Your plan for the remaining technical tasks
*   External tools you are using or plan to use

## 2 Mid-term Oral Report

You are going to attend the session scheduled on Dec 4, 5-6, to show your project progress **in person** with demos you might have already finished. It is not necessary to have completed all features by this time, but make sure to have a working demo for presentation. Notice that the oral report is an important criterion for the TAs to select orals for the final presentation.

## 3 In-class Presentation

Each team is required to create a poster and provide live demonstrations for the final presentation. The three TAs and the instructor will evaluate each poster based on the quality of the demo—considering factors like the aesthetic appeal of your rendering, the gaming experience of your interactive game, or the complexity of your simulation—as well as the overall technical presentation and the effectiveness of the on-site Q&A.

Additionally, all enrolled students will have the chance to vote for their favorite projects, with the most popular projects earning up to 3 bonus points. We will also award 2 bonus points for oral presentations. The TAs will select approximately 10 projects for oral presentations based on the mid-term oral report. In addition to the required poster presentation, these selected projects will have the opportunity to be presented orally in class on December 24.

## 4 Final Written Report

Each individual must submit a final report detailing their personal contributions using the SIGGRAPH template. The report should be at least three pages long and include the following sections: Introduction, Method, Results, Discussion, Personal Contribution Statement, and References. Each report must be written independently—**no shared writing**, even within teams.

In addition, you must submit the source code for the project (a GitHub link is recommended). Whether working individually or in a team, proper commit management is advised to help TAs evaluate your workload. In the final written report, make sure clearly indicate which part of the code you contributed.

## 5 Projects

In the following sections, the **red parts** represent the core requirements for this topic. You must **fully implement** these components; failure to do so will result in your project being considered incomplete. It is acceptable if you choose not to pursue the optional parts, but skipping the entire project is **not allowed**.

There are two possible scenarios:

*   **Individual:** You must complete the **red parts**, which will automatically earn you 5pts. Afterward, you may optionally select any combination of the **black parts** to earn up to an additional 10pts.
*   **Team:** Your team must complete the **red parts**, which will automatically earn each individual 8pts. Then, each individual may optionally select a **different** combination of the **black parts** to earn up to an additional 7pts. This means the red parts are completed collaboratively, but only one team member can receive points for each of the black parts.

Overall, we encourage everyone to complete the project as a team. Since the high difficulty of some technical points, we believe that completing 10pts individually is comparable in difficulty to completing 14pts as a team, so the total workload is similar. However, if you choose to work as a team, you are more likely to achieve higher completion and presentation quality, which can help you earn better scores and possibly bonuses during the presentation portion.

For the writing portion, the TA will evaluate the quality of each individual’s report. Deductions of 1-5 pts may be applied for issues such as improper structure, unclear methodology, missing experimental results, or insufficient length.

If you wish to implement any advanced techniques in your project, please contact the TA before the **Oral Report** to discuss feasibility and grading criteria.

### 5.1 Topic 1: Image Rendering

In this section, you are required to implement a renderer, which can be either hardware-based (GPU) or software-based (CPU). The goal is to produce a static image. Note that a path-tracing renderer is required, as rasterization-based renderers will not be accepted. Also, there is a framework that you can choose to use it or not (github.com/LazyJazzDev/Sparks). The functionalities you need to implement are as follows:

*   **Base:** Implement a path tracing algorithm that correctly handles **diffuse** and **specular** materials. (basic)
*   **Scene creation:** Build a custom scene with aesthetic considerations, using geometry that you create from scratch or find online (ensure the source is credited). (basic, tidiness and attractiveness 1pt)
*   **Acceleration structure:** Implement an acceleration structure such as BVH (Bounding Volume Hierarchy). This is not required for hardware-based renderers, as the acceleration structure is built-in in that case. (basic, Surface Area Heuristic or another advanced algorithm 2pts)
*   **Material:** Create a (non-trivial) custom material. Options include:
    *   Transmissive material (basic)
    *   Principled BSDF (2pts)
    *   Multi-layer material (2pts)
    *   Rendering of fur, hair, skin, etc. (2pts)
*   **Texture:** Create your own (non-trivial) texture with proper texture mapping. Options include:
    *   Color texture (basic)
    *   Normal map, height map, attribute map, or any functional texture mapping (1pt for each, up to 2pts)
    *   Implement an adaptive mipmap algorithm (2pts)
*   **Importance Sampling:** Use more advanced sampling algorithms for path tracing. (Importance sampling with Russian Roulette, multiple importance sampling 2pts)
*   **Volumetric Rendering:** Options include:
    *   Subsurface scattering (2pts)
    *   Homogeneous volume rendering (1pt)
    *   Inhomogeneous volume rendering (1pt)
    *   Channel-independent subsurface scattering (1pt)
    *   Volumetric emission (1pt)
    *   Volumetric alpha shadow (2pts)
*   **Special Visual Effects:** Options include:
    *   Motion blur, depth of field (basic)
    *   Alpha shadow (basic)
    *   Cartoon style rendering (2pts)
    *   Chromatic dispersion (2pts)
*   **Lighting:** Options include:
    *   Point light and area light (basic)
    *   Environment lighting with HDR, such as skybox (2pts)
*   **Anti-aliasing:** Implement an anti-aliasing algorithm (basic)
*   **Simulation-based content creation:** Up to 2pts